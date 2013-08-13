;(function(e,t,n){function i(n,s){if(!t[n]){if(!e[n]){var o=typeof require=="function"&&require;if(!s&&o)return o(n,!0);if(r)return r(n,!0);throw new Error("Cannot find module '"+n+"'")}var u=t[n]={exports:{}};e[n][0].call(u.exports,function(t){var r=e[n][1][t];return i(r?r:t)},u,u.exports)}return t[n].exports}var r=typeof require=="function"&&require;for(var s=0;s<n.length;s++)i(n[s]);return i})({1:[function(require,module,exports){

/*

Deterministic fire model based on firelib + a cellular automata 
fire growth model (FGM)

template fuction receives a three element array "dataArray" with moisture[%], wind speed [m/s]
and wind direction[º from north] 

.other variables:
  slopeMapPC   - Slope Map array
  aspectMapPC  - Aspect Map Array
  clcMapPC     - Corine land cover map 
  rowsPC       - Number of rows
  colsPC       - Number os Columnss
  heightPC     - Map height
  widthPC      - Map width

*/


//!!!ACHTUNG - Don't Fuck with the fuel model. 

module.exports = function (dataArray, rowsPC, colsPC, aspectMapPC, slopeMapPC, clcMapPC, heightPC, widthPC){

  var fireLib = require('./fireLib');
  //var fireLib = require('./slowFGM');

  var rows = rowsPC;
  var cols = colsPC;
  var MOISTUREPART = dataArray[0]/100;             //fraction
  var WINDU = dataArray[1]*196.850393701;          // [m/s] - > ft/min (2.23 m/s = 5mph)
  var WINDDIR =dataArray[2];                       //degrees clockwise from north


  var H = metersToFeet(heightPC);                      //Terrain Length
  var W = metersToFeet(widthPC);                       //Terrain Width

  var CellWd = W/rows;
  var CellHt = H/cols;

  var INF = 9999999999999;
  var smidgen = 1E-6;

  var nStencil = 16;

  var row, col, nrow, ncol;
  var cell;
  var cells = rows*cols;
  var ncell;
  var dCell;
  
             //N   NE   E  SE  S  SW   W  NW   a   b   c   d   e  f   g  h 
  var nCol = [ 0,   1,  1,  1, 0, -1, -1, -1, -1,  1, -2,  2, -2, 2, -1, 1];
  var nRow = [ -1, -1,  0,  1, 1,  1,  0, -1, -2, -2, -1, -1,  1, 1,  2, 2];
  var nDist = new Array (nStencil);
  var nAzm =  new Array (nStencil);

  var timeNext = 0;
  var timeNow = 0;
  var ignNcell;
  var ignTime;

  //create maps
  var ignMap            = new Array (rows*cols);
  var ignMapNew         = new Array (rows*cols);    //Used in iterative (Fast) FGM
  var rosMap            = new Array (rows*cols);
  var rosMaxMap         = new Array (rows*cols);
  var ros0Map           = new Array (rows*cols);
  var rxIntensityMap    = new Array (rows*cols);
  var moistMap          = new Array (rows*cols); 
  var windUMap          = new Array (rows*cols); 
  var windDirMap        = new Array (rows*cols); 
  var slopeMap          = new Array (rows*cols);
  var aspectMap         = new Array (rows*cols);
  var phiEffWindMap     = new Array (rows*cols);
  var eccentricityMap   = new Array (rows*cols);
  var azimuthMaxMap     = new Array (rows*cols);

  //Read file properties, build fuelProps object
  var fuelProps = createFuelPropsNFFL1();


  initMaps();

  FGM();

  for (cell = 0; cell < rows*cols; cell++)
    ignMap[cell] = parseFloat(ignMap[cell].toFixed(2));

  return JSON.stringify(ignMap);

  function FGM(){

    //Compute distance and Azimuth of neighbour
    //in a outward propagation configuration
    calcDistAzm();


    while (timeNext < INF){
      timeNow = timeNext;
      timeNext = INF;

      for ( row = 0; row < rows; row++){
        for ( col = 0; col < cols; col++){
          cell = col + cols*row;
          
          //If the cell burns only in the future, skips and update timeNext if necessary
          //finds the minimum timeNext from the cells ignition times
          if ( ignMap[cell] > timeNow && timeNext > ignMap[cell] ){

            timeNext = ignMap[cell];
            continue;
          } 
          if ( ignMap[cell] !== timeNow )
            continue;

          //Neighbour loop if ignMap[cell] = timeNow
          for (var n = 0; n < 16; n++){

            //neighbour index calc
            ncol = col + nCol[n];
            nrow = row + nRow[n];
            ncell = ncol + nrow*cols;

            //Check if neighbour is inbound
            if ( !(nrow >= 0 && nrow < rows && ncol >= 0 && ncol < cols))
              continue;


            var ignNcell = ignMap[ncell];

            // if cell is unburned, compute propagation time
            if ( !(ignNcell > timeNow && rosMaxMap[cell] >= smidgen ))
              continue;

            ros = fireLib.spreadAnyAzimuth(cell, nAzm[n], phiEffWindMap, azimuthMaxMap, rosMaxMap, 
                            eccentricityMap, ros0Map );

            ignTime = timeNow + nDist[n] / ros;

            //Update ignition time
            if(ignTime < ignNcell)
              ignMap[ncell] = ignTime;

            //Update timeNext
            if( ignTime < timeNext )
              timeNext = ignTime;
          }
        }
      }
    }

    function calcDistAzm(){
      for ( n = 0; n<nStencil; n++ ){
          nDist[n] = Math.sqrt ( nCol[n] * CellWd * nCol[n] * CellWd + nRow[n] * CellHt * nRow[n] * CellHt );

          if (n < 8)
            nAzm[n] = n * 45.0;
          else
          {

            nAzm[n] = Math.atan( (nCol[n] * CellWd) / (nRow[n] * CellHt) );

            if ( nCol[n] > 0  && nRow[n] < 0) //1st quadrant 
              nAzm[n] = RadToDeg(  Math.abs( nAzm[n] ));

            if ( nCol[n] > 0  && nRow[n] > 0) //2st quadrant 
              nAzm[n] = 180.0 - RadToDeg( nAzm[n] ) ;

            if ( nCol[n] < 0  && nRow[n] > 0) //3st quadrant 
              nAzm[n] = RadToDeg( Math.abs( nAzm[n] ) )+ 180.0;

            if ( nCol[n] < 0  && nRow[n] < 0) //4st quadrant 
              nAzm[n] = 360.0 - RadToDeg( Math.abs( nAzm[n] ));
          }
      }
    }

  }

  function time(func){
    var start = Date.now();
    func();
    var end = Date.now();
    return end - start;
  }

function createFuelPropsNFFL1(){
    var array;
    var fuelObj = {};

    fuelObj.Fuel_AreaWtg = 1.00000e+00;
    fuelObj.Fuel_LifeRxFactor =1.52283e+03;
    fuelObj.Fuel_PropFlux = 5.77522e-02;
    fuelObj.Fuel_Mext = 1.20000e-01;
    fuelObj.Fuel_LifeAreaWtg = 1.00000e+00;
    fuelObj.Fuel_SigmaFactor = 9.61339e-01 ;
    fuelObj.Fuel_BulkDensity = 3.40000e-02 ;
    fuelObj.Fuel_WindB = 2.07124e+00 ;
    fuelObj.Fuel_WindK = 7.17344e-05;
    fuelObj.Fuel_SlopeK = 4.11456e+01;
    fuelObj.Fuel_WindE = 1.39403e+04;

    return fuelObj;
  }

  function createFuelPropsCustom(){
    var array;
    var fuelObj = {};

    fuelObj.Fuel_AreaWtg = 1.00000e+00;
    fuelObj.Fuel_LifeRxFactor =2.85775e+03;
    fuelObj.Fuel_PropFlux = 2.00330e+00;
    fuelObj.Fuel_Mext = 1.20000e-01;
    fuelObj.Fuel_LifeAreaWtg = 1.00000e+00;
    fuelObj.Fuel_SigmaFactor = 9.82898e-01;
    fuelObj.Fuel_BulkDensity = 1.16751e+00;
    fuelObj.Fuel_WindB = 3.23670e+00;
    fuelObj.Fuel_WindK = 5.32355e-08;
    fuelObj.Fuel_SlopeK = 1.42426e+01;
    fuelObj.Fuel_WindE = 1.87845e+07;

    return fuelObj;
  }

  function initMaps(){

    //Init maps
    for (cell = 0; cell < cells; cell++){
      ignMap[cell]      = INF;
      moistMap[cell]    = MOISTUREPART;
      windUMap[cell]    = WINDU;
      windDirMap[cell]  = WINDDIR;
      //Aspect in firelib is N=0 and clockwise 
      //while aspect in Grass is E=0 counter-clockwise
      aspectMap[cell] = (aspectMapPC[cell] - 90 < 0) ?                            
                          aspectMapPC[cell] - 90 + 360  : aspectMapPC[cell] - 90 ; 
      aspectMap[cell] = 360 - aspectMap[cell];
      //while in Grass is percentage rise/reach.

      //Slope in firelib is a fraction
      slopeMap[cell]    = slopeMapPC[cell]/100;

    }

    for (cell = 0; cell < cells; cell++)
      ros0Map[cell] = fireLib.noWindNoSlope(cell, fuelProps, moistMap, rxIntensityMap);


    for (cell = 0; cell < cells; cell++)
      rosMaxMap[cell] = fireLib.windAndSlope(cell, fuelProps, slopeMap, ros0Map, windUMap,
                        windDirMap, aspectMap, azimuthMaxMap, eccentricityMap,
                        phiEffWindMap, rxIntensityMap);


    //Ignition point at terrain midle
    ignMap[Math.floor(cols/2) + Math.floor(rows/2)*cols] = 0;

    //the clc maps are used to decide if a cell is burnable or not
    //Every clc value equal to 32X or 31X is considered to be a custom fuel model
    //otherwise Ros0 and RosMax are zero
    for (var n = 0; n<clcMapPC.length; n++){
      if (  !(/\b32\d\b/.test(clcMapPC[n]) || /\b31\d\b/.test(clcMapPC[n]))  ){
        ros0Map[n] = 0;
        rosMaxMap[n] = 0;
      }
    }
  }

  function feetToMeters(x){
    x *= 0.3048;
    return x;
  }

  function metersToFeet(x){
    x *= 3.2808399;
    return x;
  }

  function DegToRad(x) {
    x *= 0.017453293;
    return x;
  }

  function RadToDeg(x) {
    x *= 57.29577951;
    return x;
  }
};
},{"./fireLib":2}],2:[function(require,module,exports){
/*

  firelib porting to javascript

  This should be really moved to it's own module

*/

var smidgen = 1E-6;
var M_PI = 3.141592653589793;
var INF = 9999999999999;


function noWindNoSlope(idx, fuelProps, moistMap, rxIntensityMap){

  var AreaWtg;
  var ratio; 
  var rxIntensity;
  var Spread0Idx;


  AreaWtg = fuelProps.Fuel_AreaWtg;

  ratio = AreaWtg*moistMap[idx]/fuelProps.Fuel_Mext;

  rxIntensity =  fuelProps.Fuel_LifeRxFactor*
  (1-2.59*ratio + 5.11*ratio*ratio - 3.52*ratio*ratio*ratio); //EtaM
  
  rxIntensityMap[idx] = rxIntensity;

  Spread0Idx = fuelProps.Fuel_PropFlux*rxIntensity /
                      ((250.0 + 1116.0*moistMap[idx])*AreaWtg*    //Qig - Heat of pre Ignition
                       fuelProps.Fuel_LifeAreaWtg*
                       fuelProps.Fuel_SigmaFactor*
                       fuelProps.Fuel_BulkDensity);

  return Spread0Idx;
}

function windAndSlope(idx, fuelProps, slopeMap, ros0Map, windUMap, windDirMap, aspectMap,
                        azimuthMaxMap, eccentricityMap, phiEffWindMap, rxIntensityMap )
{

  var windB, windK,  phiSlope, phiWind, phiEw, upSlope, spreadMax, spreadMaxIdx;
  var spread0Idx;
  var slope;
  var effectiveWind;
  var maxWind;
  var lwRatio;
  var split;
  var x;
  var y;
  var Rv;
  var a;

  slope  = slopeMap[idx];
  spread0Idx = ros0Map[idx];

  windB = fuelProps.Fuel_WindB;
  windK = fuelProps.Fuel_WindK;

  phiSlope = fuelProps.Fuel_SlopeK*slope *slope;
  phiWind  = fuelProps.Fuel_WindK*Math.pow(windUMap[idx],windB);

  //PhiWind tem um teste < smidgen em relacao a velocidade do vento WindUMap
  phiEw = phiSlope + phiWind;

  if((upSlope = aspectMap[idx]) >= 180.0)
    upSlope = upSlope - 180;
  else
    upSlope = upSlope + 180;


  //Situation 1 No fire Spread or reaction Intensity
  if(spread0Idx < smidgen) {
    spreadMaxIdx          = 0;
    eccentricityMap[idx]  = 0;
    azimuthMaxMap[idx]    = 0;
    phiEffWindMap[idx]    = phiEw;
  }

  //Situation 2 No Wind and No Slope
  else if (phiEw < smidgen) {
    phiEffWindMap[idx]   = 0;
    spreadMaxIdx         = spread0Idx;
    eccentricityMap[idx] = 0;
    azimuthMaxMap[idx]   = 0;
  }

  //Situation 3 Wind with No Slope
  else if (slope  < smidgen) {
    effectiveWind  = windUMap[idx];
    azimuthMaxMap[idx] = windDirMap[idx];
    
    maxWind = 0.9*rxIntensityMap[idx];
    spread0Idx = ros0Map[idx];
    if(effectiveWind  >  maxWind ) {

      phiEw = windK*Math.pow(maxWind , windB);
      effectiveWind  = maxWind ;
    }

    spreadMaxIdx = spread0Idx *(1 + phiEw);
    
    if(effectiveWind  >  smidgen) {

      lwRatio  = 1.0 + 0.002840909 * effectiveWind ;
      if (lwRatio  > 1.00001)
        eccentricityMap[idx] = Math.sqrt(lwRatio *lwRatio  - 1)/lwRatio ;
    }

    phiEffWindMap[idx]  = phiEw;
  }

  //Situation 4 and 5 - slope with no wind and wind blows upSlope
  else if(windUMap[idx] < smidgen || equal(upSlope, windDirMap[idx])) {

    azimuthMaxMap[idx] = upSlope;
    effectiveWind  = Math.pow(phiEw*fuelProps.Fuel_WindE, 1/windB);
    maxWind  = 0.9*rxIntensityMap[idx];

    if(effectiveWind  >  maxWind ) {

      phiEw = windK*Math.pow(maxWind , windB);
      effectiveWind  = maxWind ;
    }

    if(effectiveWind  >  smidgen) {

      lwRatio  = 1.0 + 0.002840909 * effectiveWind;
      if (lwRatio  > 1.000001)
        eccentricityMap[idx] = Math.sqrt(lwRatio *lwRatio  - 1)/lwRatio ;
    }

    spreadMaxIdx = spread0Idx *(1 + phiEw);
  }
  //Situation 6 - Wind Blows cross Slope
  else {

    split  = windDirMap[idx];
    if (upSlope <= split )
      split  = split  - upSlope;
    else
      split  = 360.0 - upSlope + split ;

    split  = DegToRad(split );
    x   = spread0Idx *(phiSlope + phiWind*Math.cos(split ));
    y   = spread0Idx *(phiWind*Math.sin(split ));
    Rv  = Math.sqrt(x *x  + y *y );

    spreadMax = spread0Idx  + Rv ;
    phiEw = spreadMax / spread0Idx  - 1;
    a  = Math.asin(Math.abs(y ) / Rv );
    if(x  >= 0.0)
      a  = (y  >= 0.0) ? a           : M_PI + M_PI - a ;
    else
      a  = (y  >= 0.0) ? (M_PI - a ) : M_PI + a ;
    
    split  = RadToDeg(a );
    if (upSlope + split  > 306.0)
      azimuthMaxMap[idx] = upSlope + split  - 360.0;
    else
      azimuthMaxMap[idx] = upSlope + split ;

    effectiveWind  = Math.pow(phiEw*fuelProps.Fuel_WindE, 1/windB);
    
    //Do effective wind only if phiEw > smidgen
    if(phiEw > smidgen) {

      maxWind  = 0.9*rxIntensityMap[idx];
      if(effectiveWind  >  maxWind ) {
        phiEw = windK*Math.pow(maxWind , windB);
        effectiveWind  = maxWind ;
        spreadMax = spread0Idx *(1 + phiEw);
      }
    }

    if(effectiveWind  >  smidgen) {

      lwRatio  = 1.0 + 0.002840909 * effectiveWind ;
      if (lwRatio  > 1.00001)
        eccentricityMap[idx] = Math.sqrt(lwRatio *lwRatio  - 1)/lwRatio ;
    }

    spreadMaxIdx = spreadMax;
    phiEffWindMap[idx] = phiEw;

  }

  return ( spreadMaxIdx );
}

function spreadAnyAzimuth(idx, azimuth, phiEffWindMap, azimuthMaxMap, rosMaxMap,
                            eccentricityMap, ros0Map )
{

  var spreadAny;
  var eccentricity;

  if (phiEffWindMap[idx] < smidgen && azimuthMaxMap[idx] === azimuth)
  
    spreadAny = rosMaxMap[idx];
  
  else {
  
    if ((dir = Math.abs(azimuthMaxMap[idx] - azimuth)) > 180)
      dir = 360.0 - dir;

    dir = DegToRad(dir);

    eccentricity = eccentricityMap[idx];
    spreadAny = rosMaxMap[idx]*(1 - eccentricity)/(1 - eccentricity*Math.cos(dir));

    if (spreadAny > INF)
      spreadAny = ros0Map[idx];

  }

  return spreadAny;
}


function DegToRad(x) {
  x *= 0.017453293;
  return x;
}

function RadToDeg(x) {
  x *= 57.29577951;
  return x;
}

function equal(x,y){
  if ( Math.abs(x-y)<smidgen )
    return true;
  else
    return false;
}

exports.windAndSlope = windAndSlope;
exports.spreadAnyAzimuth = spreadAnyAzimuth;
exports.noWindNoSlope = noWindNoSlope;

},{}]},{},[1])
//@ sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoiZ2VuZXJhdGVkLmpzIiwic291cmNlcyI6WyIvaG9tZS9mc291c2Evc3JjL2NycC9lbWJlcnMvZW5naW5lL3NyYy9wcm9ncmFtLmpzIiwiL2hvbWUvZnNvdXNhL3NyYy9jcnAvZW1iZXJzL2VuZ2luZS9zcmMvZmlyZUxpYi5qcyJdLCJuYW1lcyI6W10sIm1hcHBpbmdzIjoiO0FBQUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNqU0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBIiwic291cmNlc0NvbnRlbnQiOlsiXG4vKlxuXG5EZXRlcm1pbmlzdGljIGZpcmUgbW9kZWwgYmFzZWQgb24gZmlyZWxpYiArIGEgY2VsbHVsYXIgYXV0b21hdGEgXG5maXJlIGdyb3d0aCBtb2RlbCAoRkdNKVxuXG50ZW1wbGF0ZSBmdWN0aW9uIHJlY2VpdmVzIGEgdGhyZWUgZWxlbWVudCBhcnJheSBcImRhdGFBcnJheVwiIHdpdGggbW9pc3R1cmVbJV0sIHdpbmQgc3BlZWQgW20vc11cbmFuZCB3aW5kIGRpcmVjdGlvblvCuiBmcm9tIG5vcnRoXSBcblxuLm90aGVyIHZhcmlhYmxlczpcbiAgc2xvcGVNYXBQQyAgIC0gU2xvcGUgTWFwIGFycmF5XG4gIGFzcGVjdE1hcFBDICAtIEFzcGVjdCBNYXAgQXJyYXlcbiAgY2xjTWFwUEMgICAgIC0gQ29yaW5lIGxhbmQgY292ZXIgbWFwIFxuICByb3dzUEMgICAgICAgLSBOdW1iZXIgb2Ygcm93c1xuICBjb2xzUEMgICAgICAgLSBOdW1iZXIgb3MgQ29sdW1uc3NcbiAgaGVpZ2h0UEMgICAgIC0gTWFwIGhlaWdodFxuICB3aWR0aFBDICAgICAgLSBNYXAgd2lkdGhcblxuKi9cblxuXG4vLyEhIUFDSFRVTkcgLSBEb24ndCBGdWNrIHdpdGggdGhlIGZ1ZWwgbW9kZWwuIFxuXG5tb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uIChkYXRhQXJyYXksIHJvd3NQQywgY29sc1BDLCBhc3BlY3RNYXBQQywgc2xvcGVNYXBQQywgY2xjTWFwUEMsIGhlaWdodFBDLCB3aWR0aFBDKXtcblxuICB2YXIgZmlyZUxpYiA9IHJlcXVpcmUoJy4vZmlyZUxpYicpO1xuICAvL3ZhciBmaXJlTGliID0gcmVxdWlyZSgnLi9zbG93RkdNJyk7XG5cbiAgdmFyIHJvd3MgPSByb3dzUEM7XG4gIHZhciBjb2xzID0gY29sc1BDO1xuICB2YXIgTU9JU1RVUkVQQVJUID0gZGF0YUFycmF5WzBdLzEwMDsgICAgICAgICAgICAgLy9mcmFjdGlvblxuICB2YXIgV0lORFUgPSBkYXRhQXJyYXlbMV0qMTk2Ljg1MDM5MzcwMTsgICAgICAgICAgLy8gW20vc10gLSA+IGZ0L21pbiAoMi4yMyBtL3MgPSA1bXBoKVxuICB2YXIgV0lORERJUiA9ZGF0YUFycmF5WzJdOyAgICAgICAgICAgICAgICAgICAgICAgLy9kZWdyZWVzIGNsb2Nrd2lzZSBmcm9tIG5vcnRoXG5cblxuICB2YXIgSCA9IG1ldGVyc1RvRmVldChoZWlnaHRQQyk7ICAgICAgICAgICAgICAgICAgICAgIC8vVGVycmFpbiBMZW5ndGhcbiAgdmFyIFcgPSBtZXRlcnNUb0ZlZXQod2lkdGhQQyk7ICAgICAgICAgICAgICAgICAgICAgICAvL1RlcnJhaW4gV2lkdGhcblxuICB2YXIgQ2VsbFdkID0gVy9yb3dzO1xuICB2YXIgQ2VsbEh0ID0gSC9jb2xzO1xuXG4gIHZhciBJTkYgPSA5OTk5OTk5OTk5OTk5O1xuICB2YXIgc21pZGdlbiA9IDFFLTY7XG5cbiAgdmFyIG5TdGVuY2lsID0gMTY7XG5cbiAgdmFyIHJvdywgY29sLCBucm93LCBuY29sO1xuICB2YXIgY2VsbDtcbiAgdmFyIGNlbGxzID0gcm93cypjb2xzO1xuICB2YXIgbmNlbGw7XG4gIHZhciBkQ2VsbDtcbiAgXG4gICAgICAgICAgICAgLy9OICAgTkUgICBFICBTRSAgUyAgU1cgICBXICBOVyAgIGEgICBiICAgYyAgIGQgICBlICBmICAgZyAgaCBcbiAgdmFyIG5Db2wgPSBbIDAsICAgMSwgIDEsICAxLCAwLCAtMSwgLTEsIC0xLCAtMSwgIDEsIC0yLCAgMiwgLTIsIDIsIC0xLCAxXTtcbiAgdmFyIG5Sb3cgPSBbIC0xLCAtMSwgIDAsICAxLCAxLCAgMSwgIDAsIC0xLCAtMiwgLTIsIC0xLCAtMSwgIDEsIDEsICAyLCAyXTtcbiAgdmFyIG5EaXN0ID0gbmV3IEFycmF5IChuU3RlbmNpbCk7XG4gIHZhciBuQXptID0gIG5ldyBBcnJheSAoblN0ZW5jaWwpO1xuXG4gIHZhciB0aW1lTmV4dCA9IDA7XG4gIHZhciB0aW1lTm93ID0gMDtcbiAgdmFyIGlnbk5jZWxsO1xuICB2YXIgaWduVGltZTtcblxuICAvL2NyZWF0ZSBtYXBzXG4gIHZhciBpZ25NYXAgICAgICAgICAgICA9IG5ldyBBcnJheSAocm93cypjb2xzKTtcbiAgdmFyIGlnbk1hcE5ldyAgICAgICAgID0gbmV3IEFycmF5IChyb3dzKmNvbHMpOyAgICAvL1VzZWQgaW4gaXRlcmF0aXZlIChGYXN0KSBGR01cbiAgdmFyIHJvc01hcCAgICAgICAgICAgID0gbmV3IEFycmF5IChyb3dzKmNvbHMpO1xuICB2YXIgcm9zTWF4TWFwICAgICAgICAgPSBuZXcgQXJyYXkgKHJvd3MqY29scyk7XG4gIHZhciByb3MwTWFwICAgICAgICAgICA9IG5ldyBBcnJheSAocm93cypjb2xzKTtcbiAgdmFyIHJ4SW50ZW5zaXR5TWFwICAgID0gbmV3IEFycmF5IChyb3dzKmNvbHMpO1xuICB2YXIgbW9pc3RNYXAgICAgICAgICAgPSBuZXcgQXJyYXkgKHJvd3MqY29scyk7IFxuICB2YXIgd2luZFVNYXAgICAgICAgICAgPSBuZXcgQXJyYXkgKHJvd3MqY29scyk7IFxuICB2YXIgd2luZERpck1hcCAgICAgICAgPSBuZXcgQXJyYXkgKHJvd3MqY29scyk7IFxuICB2YXIgc2xvcGVNYXAgICAgICAgICAgPSBuZXcgQXJyYXkgKHJvd3MqY29scyk7XG4gIHZhciBhc3BlY3RNYXAgICAgICAgICA9IG5ldyBBcnJheSAocm93cypjb2xzKTtcbiAgdmFyIHBoaUVmZldpbmRNYXAgICAgID0gbmV3IEFycmF5IChyb3dzKmNvbHMpO1xuICB2YXIgZWNjZW50cmljaXR5TWFwICAgPSBuZXcgQXJyYXkgKHJvd3MqY29scyk7XG4gIHZhciBhemltdXRoTWF4TWFwICAgICA9IG5ldyBBcnJheSAocm93cypjb2xzKTtcblxuICAvL1JlYWQgZmlsZSBwcm9wZXJ0aWVzLCBidWlsZCBmdWVsUHJvcHMgb2JqZWN0XG4gIHZhciBmdWVsUHJvcHMgPSBjcmVhdGVGdWVsUHJvcHNORkZMMSgpO1xuXG5cbiAgaW5pdE1hcHMoKTtcblxuICBGR00oKTtcblxuICBmb3IgKGNlbGwgPSAwOyBjZWxsIDwgcm93cypjb2xzOyBjZWxsKyspXG4gICAgaWduTWFwW2NlbGxdID0gcGFyc2VGbG9hdChpZ25NYXBbY2VsbF0udG9GaXhlZCgyKSk7XG5cbiAgcmV0dXJuIEpTT04uc3RyaW5naWZ5KGlnbk1hcCk7XG5cbiAgZnVuY3Rpb24gRkdNKCl7XG5cbiAgICAvL0NvbXB1dGUgZGlzdGFuY2UgYW5kIEF6aW11dGggb2YgbmVpZ2hib3VyXG4gICAgLy9pbiBhIG91dHdhcmQgcHJvcGFnYXRpb24gY29uZmlndXJhdGlvblxuICAgIGNhbGNEaXN0QXptKCk7XG5cblxuICAgIHdoaWxlICh0aW1lTmV4dCA8IElORil7XG4gICAgICB0aW1lTm93ID0gdGltZU5leHQ7XG4gICAgICB0aW1lTmV4dCA9IElORjtcblxuICAgICAgZm9yICggcm93ID0gMDsgcm93IDwgcm93czsgcm93Kyspe1xuICAgICAgICBmb3IgKCBjb2wgPSAwOyBjb2wgPCBjb2xzOyBjb2wrKyl7XG4gICAgICAgICAgY2VsbCA9IGNvbCArIGNvbHMqcm93O1xuICAgICAgICAgIFxuICAgICAgICAgIC8vSWYgdGhlIGNlbGwgYnVybnMgb25seSBpbiB0aGUgZnV0dXJlLCBza2lwcyBhbmQgdXBkYXRlIHRpbWVOZXh0IGlmIG5lY2Vzc2FyeVxuICAgICAgICAgIC8vZmluZHMgdGhlIG1pbmltdW0gdGltZU5leHQgZnJvbSB0aGUgY2VsbHMgaWduaXRpb24gdGltZXNcbiAgICAgICAgICBpZiAoIGlnbk1hcFtjZWxsXSA+IHRpbWVOb3cgJiYgdGltZU5leHQgPiBpZ25NYXBbY2VsbF0gKXtcblxuICAgICAgICAgICAgdGltZU5leHQgPSBpZ25NYXBbY2VsbF07XG4gICAgICAgICAgICBjb250aW51ZTtcbiAgICAgICAgICB9IFxuICAgICAgICAgIGlmICggaWduTWFwW2NlbGxdICE9PSB0aW1lTm93IClcbiAgICAgICAgICAgIGNvbnRpbnVlO1xuXG4gICAgICAgICAgLy9OZWlnaGJvdXIgbG9vcCBpZiBpZ25NYXBbY2VsbF0gPSB0aW1lTm93XG4gICAgICAgICAgZm9yICh2YXIgbiA9IDA7IG4gPCAxNjsgbisrKXtcblxuICAgICAgICAgICAgLy9uZWlnaGJvdXIgaW5kZXggY2FsY1xuICAgICAgICAgICAgbmNvbCA9IGNvbCArIG5Db2xbbl07XG4gICAgICAgICAgICBucm93ID0gcm93ICsgblJvd1tuXTtcbiAgICAgICAgICAgIG5jZWxsID0gbmNvbCArIG5yb3cqY29scztcblxuICAgICAgICAgICAgLy9DaGVjayBpZiBuZWlnaGJvdXIgaXMgaW5ib3VuZFxuICAgICAgICAgICAgaWYgKCAhKG5yb3cgPj0gMCAmJiBucm93IDwgcm93cyAmJiBuY29sID49IDAgJiYgbmNvbCA8IGNvbHMpKVxuICAgICAgICAgICAgICBjb250aW51ZTtcblxuXG4gICAgICAgICAgICB2YXIgaWduTmNlbGwgPSBpZ25NYXBbbmNlbGxdO1xuXG4gICAgICAgICAgICAvLyBpZiBjZWxsIGlzIHVuYnVybmVkLCBjb21wdXRlIHByb3BhZ2F0aW9uIHRpbWVcbiAgICAgICAgICAgIGlmICggIShpZ25OY2VsbCA+IHRpbWVOb3cgJiYgcm9zTWF4TWFwW2NlbGxdID49IHNtaWRnZW4gKSlcbiAgICAgICAgICAgICAgY29udGludWU7XG5cbiAgICAgICAgICAgIHJvcyA9IGZpcmVMaWIuc3ByZWFkQW55QXppbXV0aChjZWxsLCBuQXptW25dLCBwaGlFZmZXaW5kTWFwLCBhemltdXRoTWF4TWFwLCByb3NNYXhNYXAsIFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGVjY2VudHJpY2l0eU1hcCwgcm9zME1hcCApO1xuXG4gICAgICAgICAgICBpZ25UaW1lID0gdGltZU5vdyArIG5EaXN0W25dIC8gcm9zO1xuXG4gICAgICAgICAgICAvL1VwZGF0ZSBpZ25pdGlvbiB0aW1lXG4gICAgICAgICAgICBpZihpZ25UaW1lIDwgaWduTmNlbGwpXG4gICAgICAgICAgICAgIGlnbk1hcFtuY2VsbF0gPSBpZ25UaW1lO1xuXG4gICAgICAgICAgICAvL1VwZGF0ZSB0aW1lTmV4dFxuICAgICAgICAgICAgaWYoIGlnblRpbWUgPCB0aW1lTmV4dCApXG4gICAgICAgICAgICAgIHRpbWVOZXh0ID0gaWduVGltZTtcbiAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgIH1cbiAgICB9XG5cbiAgICBmdW5jdGlvbiBjYWxjRGlzdEF6bSgpe1xuICAgICAgZm9yICggbiA9IDA7IG48blN0ZW5jaWw7IG4rKyApe1xuICAgICAgICAgIG5EaXN0W25dID0gTWF0aC5zcXJ0ICggbkNvbFtuXSAqIENlbGxXZCAqIG5Db2xbbl0gKiBDZWxsV2QgKyBuUm93W25dICogQ2VsbEh0ICogblJvd1tuXSAqIENlbGxIdCApO1xuXG4gICAgICAgICAgaWYgKG4gPCA4KVxuICAgICAgICAgICAgbkF6bVtuXSA9IG4gKiA0NS4wO1xuICAgICAgICAgIGVsc2VcbiAgICAgICAgICB7XG5cbiAgICAgICAgICAgIG5Bem1bbl0gPSBNYXRoLmF0YW4oIChuQ29sW25dICogQ2VsbFdkKSAvIChuUm93W25dICogQ2VsbEh0KSApO1xuXG4gICAgICAgICAgICBpZiAoIG5Db2xbbl0gPiAwICAmJiBuUm93W25dIDwgMCkgLy8xc3QgcXVhZHJhbnQgXG4gICAgICAgICAgICAgIG5Bem1bbl0gPSBSYWRUb0RlZyggIE1hdGguYWJzKCBuQXptW25dICkpO1xuXG4gICAgICAgICAgICBpZiAoIG5Db2xbbl0gPiAwICAmJiBuUm93W25dID4gMCkgLy8yc3QgcXVhZHJhbnQgXG4gICAgICAgICAgICAgIG5Bem1bbl0gPSAxODAuMCAtIFJhZFRvRGVnKCBuQXptW25dICkgO1xuXG4gICAgICAgICAgICBpZiAoIG5Db2xbbl0gPCAwICAmJiBuUm93W25dID4gMCkgLy8zc3QgcXVhZHJhbnQgXG4gICAgICAgICAgICAgIG5Bem1bbl0gPSBSYWRUb0RlZyggTWF0aC5hYnMoIG5Bem1bbl0gKSApKyAxODAuMDtcblxuICAgICAgICAgICAgaWYgKCBuQ29sW25dIDwgMCAgJiYgblJvd1tuXSA8IDApIC8vNHN0IHF1YWRyYW50IFxuICAgICAgICAgICAgICBuQXptW25dID0gMzYwLjAgLSBSYWRUb0RlZyggTWF0aC5hYnMoIG5Bem1bbl0gKSk7XG4gICAgICAgICAgfVxuICAgICAgfVxuICAgIH1cblxuICB9XG5cbiAgZnVuY3Rpb24gdGltZShmdW5jKXtcbiAgICB2YXIgc3RhcnQgPSBEYXRlLm5vdygpO1xuICAgIGZ1bmMoKTtcbiAgICB2YXIgZW5kID0gRGF0ZS5ub3coKTtcbiAgICByZXR1cm4gZW5kIC0gc3RhcnQ7XG4gIH1cblxuZnVuY3Rpb24gY3JlYXRlRnVlbFByb3BzTkZGTDEoKXtcbiAgICB2YXIgYXJyYXk7XG4gICAgdmFyIGZ1ZWxPYmogPSB7fTtcblxuICAgIGZ1ZWxPYmouRnVlbF9BcmVhV3RnID0gMS4wMDAwMGUrMDA7XG4gICAgZnVlbE9iai5GdWVsX0xpZmVSeEZhY3RvciA9MS41MjI4M2UrMDM7XG4gICAgZnVlbE9iai5GdWVsX1Byb3BGbHV4ID0gNS43NzUyMmUtMDI7XG4gICAgZnVlbE9iai5GdWVsX01leHQgPSAxLjIwMDAwZS0wMTtcbiAgICBmdWVsT2JqLkZ1ZWxfTGlmZUFyZWFXdGcgPSAxLjAwMDAwZSswMDtcbiAgICBmdWVsT2JqLkZ1ZWxfU2lnbWFGYWN0b3IgPSA5LjYxMzM5ZS0wMSA7XG4gICAgZnVlbE9iai5GdWVsX0J1bGtEZW5zaXR5ID0gMy40MDAwMGUtMDIgO1xuICAgIGZ1ZWxPYmouRnVlbF9XaW5kQiA9IDIuMDcxMjRlKzAwIDtcbiAgICBmdWVsT2JqLkZ1ZWxfV2luZEsgPSA3LjE3MzQ0ZS0wNTtcbiAgICBmdWVsT2JqLkZ1ZWxfU2xvcGVLID0gNC4xMTQ1NmUrMDE7XG4gICAgZnVlbE9iai5GdWVsX1dpbmRFID0gMS4zOTQwM2UrMDQ7XG5cbiAgICByZXR1cm4gZnVlbE9iajtcbiAgfVxuXG4gIGZ1bmN0aW9uIGNyZWF0ZUZ1ZWxQcm9wc0N1c3RvbSgpe1xuICAgIHZhciBhcnJheTtcbiAgICB2YXIgZnVlbE9iaiA9IHt9O1xuXG4gICAgZnVlbE9iai5GdWVsX0FyZWFXdGcgPSAxLjAwMDAwZSswMDtcbiAgICBmdWVsT2JqLkZ1ZWxfTGlmZVJ4RmFjdG9yID0yLjg1Nzc1ZSswMztcbiAgICBmdWVsT2JqLkZ1ZWxfUHJvcEZsdXggPSAyLjAwMzMwZSswMDtcbiAgICBmdWVsT2JqLkZ1ZWxfTWV4dCA9IDEuMjAwMDBlLTAxO1xuICAgIGZ1ZWxPYmouRnVlbF9MaWZlQXJlYVd0ZyA9IDEuMDAwMDBlKzAwO1xuICAgIGZ1ZWxPYmouRnVlbF9TaWdtYUZhY3RvciA9IDkuODI4OThlLTAxO1xuICAgIGZ1ZWxPYmouRnVlbF9CdWxrRGVuc2l0eSA9IDEuMTY3NTFlKzAwO1xuICAgIGZ1ZWxPYmouRnVlbF9XaW5kQiA9IDMuMjM2NzBlKzAwO1xuICAgIGZ1ZWxPYmouRnVlbF9XaW5kSyA9IDUuMzIzNTVlLTA4O1xuICAgIGZ1ZWxPYmouRnVlbF9TbG9wZUsgPSAxLjQyNDI2ZSswMTtcbiAgICBmdWVsT2JqLkZ1ZWxfV2luZEUgPSAxLjg3ODQ1ZSswNztcblxuICAgIHJldHVybiBmdWVsT2JqO1xuICB9XG5cbiAgZnVuY3Rpb24gaW5pdE1hcHMoKXtcblxuICAgIC8vSW5pdCBtYXBzXG4gICAgZm9yIChjZWxsID0gMDsgY2VsbCA8IGNlbGxzOyBjZWxsKyspe1xuICAgICAgaWduTWFwW2NlbGxdICAgICAgPSBJTkY7XG4gICAgICBtb2lzdE1hcFtjZWxsXSAgICA9IE1PSVNUVVJFUEFSVDtcbiAgICAgIHdpbmRVTWFwW2NlbGxdICAgID0gV0lORFU7XG4gICAgICB3aW5kRGlyTWFwW2NlbGxdICA9IFdJTkRESVI7XG4gICAgICAvL0FzcGVjdCBpbiBmaXJlbGliIGlzIE49MCBhbmQgY2xvY2t3aXNlIFxuICAgICAgLy93aGlsZSBhc3BlY3QgaW4gR3Jhc3MgaXMgRT0wIGNvdW50ZXItY2xvY2t3aXNlXG4gICAgICBhc3BlY3RNYXBbY2VsbF0gPSAoYXNwZWN0TWFwUENbY2VsbF0gLSA5MCA8IDApID8gICAgICAgICAgICAgICAgICAgICAgICAgICAgXG4gICAgICAgICAgICAgICAgICAgICAgICAgIGFzcGVjdE1hcFBDW2NlbGxdIC0gOTAgKyAzNjAgIDogYXNwZWN0TWFwUENbY2VsbF0gLSA5MCA7IFxuICAgICAgYXNwZWN0TWFwW2NlbGxdID0gMzYwIC0gYXNwZWN0TWFwW2NlbGxdO1xuICAgICAgLy93aGlsZSBpbiBHcmFzcyBpcyBwZXJjZW50YWdlIHJpc2UvcmVhY2guXG5cbiAgICAgIC8vU2xvcGUgaW4gZmlyZWxpYiBpcyBhIGZyYWN0aW9uXG4gICAgICBzbG9wZU1hcFtjZWxsXSAgICA9IHNsb3BlTWFwUENbY2VsbF0vMTAwO1xuXG4gICAgfVxuXG4gICAgZm9yIChjZWxsID0gMDsgY2VsbCA8IGNlbGxzOyBjZWxsKyspXG4gICAgICByb3MwTWFwW2NlbGxdID0gZmlyZUxpYi5ub1dpbmROb1Nsb3BlKGNlbGwsIGZ1ZWxQcm9wcywgbW9pc3RNYXAsIHJ4SW50ZW5zaXR5TWFwKTtcblxuXG4gICAgZm9yIChjZWxsID0gMDsgY2VsbCA8IGNlbGxzOyBjZWxsKyspXG4gICAgICByb3NNYXhNYXBbY2VsbF0gPSBmaXJlTGliLndpbmRBbmRTbG9wZShjZWxsLCBmdWVsUHJvcHMsIHNsb3BlTWFwLCByb3MwTWFwLCB3aW5kVU1hcCxcbiAgICAgICAgICAgICAgICAgICAgICAgIHdpbmREaXJNYXAsIGFzcGVjdE1hcCwgYXppbXV0aE1heE1hcCwgZWNjZW50cmljaXR5TWFwLFxuICAgICAgICAgICAgICAgICAgICAgICAgcGhpRWZmV2luZE1hcCwgcnhJbnRlbnNpdHlNYXApO1xuXG5cbiAgICAvL0lnbml0aW9uIHBvaW50IGF0IHRlcnJhaW4gbWlkbGVcbiAgICBpZ25NYXBbTWF0aC5mbG9vcihjb2xzLzIpICsgTWF0aC5mbG9vcihyb3dzLzIpKmNvbHNdID0gMDtcblxuICAgIC8vdGhlIGNsYyBtYXBzIGFyZSB1c2VkIHRvIGRlY2lkZSBpZiBhIGNlbGwgaXMgYnVybmFibGUgb3Igbm90XG4gICAgLy9FdmVyeSBjbGMgdmFsdWUgZXF1YWwgdG8gMzJYIG9yIDMxWCBpcyBjb25zaWRlcmVkIHRvIGJlIGEgY3VzdG9tIGZ1ZWwgbW9kZWxcbiAgICAvL290aGVyd2lzZSBSb3MwIGFuZCBSb3NNYXggYXJlIHplcm9cbiAgICBmb3IgKHZhciBuID0gMDsgbjxjbGNNYXBQQy5sZW5ndGg7IG4rKyl7XG4gICAgICBpZiAoICAhKC9cXGIzMlxcZFxcYi8udGVzdChjbGNNYXBQQ1tuXSkgfHwgL1xcYjMxXFxkXFxiLy50ZXN0KGNsY01hcFBDW25dKSkgICl7XG4gICAgICAgIHJvczBNYXBbbl0gPSAwO1xuICAgICAgICByb3NNYXhNYXBbbl0gPSAwO1xuICAgICAgfVxuICAgIH1cbiAgfVxuXG4gIGZ1bmN0aW9uIGZlZXRUb01ldGVycyh4KXtcbiAgICB4ICo9IDAuMzA0ODtcbiAgICByZXR1cm4geDtcbiAgfVxuXG4gIGZ1bmN0aW9uIG1ldGVyc1RvRmVldCh4KXtcbiAgICB4ICo9IDMuMjgwODM5OTtcbiAgICByZXR1cm4geDtcbiAgfVxuXG4gIGZ1bmN0aW9uIERlZ1RvUmFkKHgpIHtcbiAgICB4ICo9IDAuMDE3NDUzMjkzO1xuICAgIHJldHVybiB4O1xuICB9XG5cbiAgZnVuY3Rpb24gUmFkVG9EZWcoeCkge1xuICAgIHggKj0gNTcuMjk1Nzc5NTE7XG4gICAgcmV0dXJuIHg7XG4gIH1cbn07IiwiLypcblxuICBmaXJlbGliIHBvcnRpbmcgdG8gamF2YXNjcmlwdFxuXG4gIFRoaXMgc2hvdWxkIGJlIHJlYWxseSBtb3ZlZCB0byBpdCdzIG93biBtb2R1bGVcblxuKi9cblxudmFyIHNtaWRnZW4gPSAxRS02O1xudmFyIE1fUEkgPSAzLjE0MTU5MjY1MzU4OTc5MztcbnZhciBJTkYgPSA5OTk5OTk5OTk5OTk5O1xuXG5cbmZ1bmN0aW9uIG5vV2luZE5vU2xvcGUoaWR4LCBmdWVsUHJvcHMsIG1vaXN0TWFwLCByeEludGVuc2l0eU1hcCl7XG5cbiAgdmFyIEFyZWFXdGc7XG4gIHZhciByYXRpbzsgXG4gIHZhciByeEludGVuc2l0eTtcbiAgdmFyIFNwcmVhZDBJZHg7XG5cblxuICBBcmVhV3RnID0gZnVlbFByb3BzLkZ1ZWxfQXJlYVd0ZztcblxuICByYXRpbyA9IEFyZWFXdGcqbW9pc3RNYXBbaWR4XS9mdWVsUHJvcHMuRnVlbF9NZXh0O1xuXG4gIHJ4SW50ZW5zaXR5ID0gIGZ1ZWxQcm9wcy5GdWVsX0xpZmVSeEZhY3RvcipcbiAgKDEtMi41OSpyYXRpbyArIDUuMTEqcmF0aW8qcmF0aW8gLSAzLjUyKnJhdGlvKnJhdGlvKnJhdGlvKTsgLy9FdGFNXG4gIFxuICByeEludGVuc2l0eU1hcFtpZHhdID0gcnhJbnRlbnNpdHk7XG5cbiAgU3ByZWFkMElkeCA9IGZ1ZWxQcm9wcy5GdWVsX1Byb3BGbHV4KnJ4SW50ZW5zaXR5IC9cbiAgICAgICAgICAgICAgICAgICAgICAoKDI1MC4wICsgMTExNi4wKm1vaXN0TWFwW2lkeF0pKkFyZWFXdGcqICAgIC8vUWlnIC0gSGVhdCBvZiBwcmUgSWduaXRpb25cbiAgICAgICAgICAgICAgICAgICAgICAgZnVlbFByb3BzLkZ1ZWxfTGlmZUFyZWFXdGcqXG4gICAgICAgICAgICAgICAgICAgICAgIGZ1ZWxQcm9wcy5GdWVsX1NpZ21hRmFjdG9yKlxuICAgICAgICAgICAgICAgICAgICAgICBmdWVsUHJvcHMuRnVlbF9CdWxrRGVuc2l0eSk7XG5cbiAgcmV0dXJuIFNwcmVhZDBJZHg7XG59XG5cbmZ1bmN0aW9uIHdpbmRBbmRTbG9wZShpZHgsIGZ1ZWxQcm9wcywgc2xvcGVNYXAsIHJvczBNYXAsIHdpbmRVTWFwLCB3aW5kRGlyTWFwLCBhc3BlY3RNYXAsXG4gICAgICAgICAgICAgICAgICAgICAgICBhemltdXRoTWF4TWFwLCBlY2NlbnRyaWNpdHlNYXAsIHBoaUVmZldpbmRNYXAsIHJ4SW50ZW5zaXR5TWFwIClcbntcblxuICB2YXIgd2luZEIsIHdpbmRLLCAgcGhpU2xvcGUsIHBoaVdpbmQsIHBoaUV3LCB1cFNsb3BlLCBzcHJlYWRNYXgsIHNwcmVhZE1heElkeDtcbiAgdmFyIHNwcmVhZDBJZHg7XG4gIHZhciBzbG9wZTtcbiAgdmFyIGVmZmVjdGl2ZVdpbmQ7XG4gIHZhciBtYXhXaW5kO1xuICB2YXIgbHdSYXRpbztcbiAgdmFyIHNwbGl0O1xuICB2YXIgeDtcbiAgdmFyIHk7XG4gIHZhciBSdjtcbiAgdmFyIGE7XG5cbiAgc2xvcGUgID0gc2xvcGVNYXBbaWR4XTtcbiAgc3ByZWFkMElkeCA9IHJvczBNYXBbaWR4XTtcblxuICB3aW5kQiA9IGZ1ZWxQcm9wcy5GdWVsX1dpbmRCO1xuICB3aW5kSyA9IGZ1ZWxQcm9wcy5GdWVsX1dpbmRLO1xuXG4gIHBoaVNsb3BlID0gZnVlbFByb3BzLkZ1ZWxfU2xvcGVLKnNsb3BlICpzbG9wZTtcbiAgcGhpV2luZCAgPSBmdWVsUHJvcHMuRnVlbF9XaW5kSypNYXRoLnBvdyh3aW5kVU1hcFtpZHhdLHdpbmRCKTtcblxuICAvL1BoaVdpbmQgdGVtIHVtIHRlc3RlIDwgc21pZGdlbiBlbSByZWxhY2FvIGEgdmVsb2NpZGFkZSBkbyB2ZW50byBXaW5kVU1hcFxuICBwaGlFdyA9IHBoaVNsb3BlICsgcGhpV2luZDtcblxuICBpZigodXBTbG9wZSA9IGFzcGVjdE1hcFtpZHhdKSA+PSAxODAuMClcbiAgICB1cFNsb3BlID0gdXBTbG9wZSAtIDE4MDtcbiAgZWxzZVxuICAgIHVwU2xvcGUgPSB1cFNsb3BlICsgMTgwO1xuXG5cbiAgLy9TaXR1YXRpb24gMSBObyBmaXJlIFNwcmVhZCBvciByZWFjdGlvbiBJbnRlbnNpdHlcbiAgaWYoc3ByZWFkMElkeCA8IHNtaWRnZW4pIHtcbiAgICBzcHJlYWRNYXhJZHggICAgICAgICAgPSAwO1xuICAgIGVjY2VudHJpY2l0eU1hcFtpZHhdICA9IDA7XG4gICAgYXppbXV0aE1heE1hcFtpZHhdICAgID0gMDtcbiAgICBwaGlFZmZXaW5kTWFwW2lkeF0gICAgPSBwaGlFdztcbiAgfVxuXG4gIC8vU2l0dWF0aW9uIDIgTm8gV2luZCBhbmQgTm8gU2xvcGVcbiAgZWxzZSBpZiAocGhpRXcgPCBzbWlkZ2VuKSB7XG4gICAgcGhpRWZmV2luZE1hcFtpZHhdICAgPSAwO1xuICAgIHNwcmVhZE1heElkeCAgICAgICAgID0gc3ByZWFkMElkeDtcbiAgICBlY2NlbnRyaWNpdHlNYXBbaWR4XSA9IDA7XG4gICAgYXppbXV0aE1heE1hcFtpZHhdICAgPSAwO1xuICB9XG5cbiAgLy9TaXR1YXRpb24gMyBXaW5kIHdpdGggTm8gU2xvcGVcbiAgZWxzZSBpZiAoc2xvcGUgIDwgc21pZGdlbikge1xuICAgIGVmZmVjdGl2ZVdpbmQgID0gd2luZFVNYXBbaWR4XTtcbiAgICBhemltdXRoTWF4TWFwW2lkeF0gPSB3aW5kRGlyTWFwW2lkeF07XG4gICAgXG4gICAgbWF4V2luZCA9IDAuOSpyeEludGVuc2l0eU1hcFtpZHhdO1xuICAgIHNwcmVhZDBJZHggPSByb3MwTWFwW2lkeF07XG4gICAgaWYoZWZmZWN0aXZlV2luZCAgPiAgbWF4V2luZCApIHtcblxuICAgICAgcGhpRXcgPSB3aW5kSypNYXRoLnBvdyhtYXhXaW5kICwgd2luZEIpO1xuICAgICAgZWZmZWN0aXZlV2luZCAgPSBtYXhXaW5kIDtcbiAgICB9XG5cbiAgICBzcHJlYWRNYXhJZHggPSBzcHJlYWQwSWR4ICooMSArIHBoaUV3KTtcbiAgICBcbiAgICBpZihlZmZlY3RpdmVXaW5kICA+ICBzbWlkZ2VuKSB7XG5cbiAgICAgIGx3UmF0aW8gID0gMS4wICsgMC4wMDI4NDA5MDkgKiBlZmZlY3RpdmVXaW5kIDtcbiAgICAgIGlmIChsd1JhdGlvICA+IDEuMDAwMDEpXG4gICAgICAgIGVjY2VudHJpY2l0eU1hcFtpZHhdID0gTWF0aC5zcXJ0KGx3UmF0aW8gKmx3UmF0aW8gIC0gMSkvbHdSYXRpbyA7XG4gICAgfVxuXG4gICAgcGhpRWZmV2luZE1hcFtpZHhdICA9IHBoaUV3O1xuICB9XG5cbiAgLy9TaXR1YXRpb24gNCBhbmQgNSAtIHNsb3BlIHdpdGggbm8gd2luZCBhbmQgd2luZCBibG93cyB1cFNsb3BlXG4gIGVsc2UgaWYod2luZFVNYXBbaWR4XSA8IHNtaWRnZW4gfHwgZXF1YWwodXBTbG9wZSwgd2luZERpck1hcFtpZHhdKSkge1xuXG4gICAgYXppbXV0aE1heE1hcFtpZHhdID0gdXBTbG9wZTtcbiAgICBlZmZlY3RpdmVXaW5kICA9IE1hdGgucG93KHBoaUV3KmZ1ZWxQcm9wcy5GdWVsX1dpbmRFLCAxL3dpbmRCKTtcbiAgICBtYXhXaW5kICA9IDAuOSpyeEludGVuc2l0eU1hcFtpZHhdO1xuXG4gICAgaWYoZWZmZWN0aXZlV2luZCAgPiAgbWF4V2luZCApIHtcblxuICAgICAgcGhpRXcgPSB3aW5kSypNYXRoLnBvdyhtYXhXaW5kICwgd2luZEIpO1xuICAgICAgZWZmZWN0aXZlV2luZCAgPSBtYXhXaW5kIDtcbiAgICB9XG5cbiAgICBpZihlZmZlY3RpdmVXaW5kICA+ICBzbWlkZ2VuKSB7XG5cbiAgICAgIGx3UmF0aW8gID0gMS4wICsgMC4wMDI4NDA5MDkgKiBlZmZlY3RpdmVXaW5kO1xuICAgICAgaWYgKGx3UmF0aW8gID4gMS4wMDAwMDEpXG4gICAgICAgIGVjY2VudHJpY2l0eU1hcFtpZHhdID0gTWF0aC5zcXJ0KGx3UmF0aW8gKmx3UmF0aW8gIC0gMSkvbHdSYXRpbyA7XG4gICAgfVxuXG4gICAgc3ByZWFkTWF4SWR4ID0gc3ByZWFkMElkeCAqKDEgKyBwaGlFdyk7XG4gIH1cbiAgLy9TaXR1YXRpb24gNiAtIFdpbmQgQmxvd3MgY3Jvc3MgU2xvcGVcbiAgZWxzZSB7XG5cbiAgICBzcGxpdCAgPSB3aW5kRGlyTWFwW2lkeF07XG4gICAgaWYgKHVwU2xvcGUgPD0gc3BsaXQgKVxuICAgICAgc3BsaXQgID0gc3BsaXQgIC0gdXBTbG9wZTtcbiAgICBlbHNlXG4gICAgICBzcGxpdCAgPSAzNjAuMCAtIHVwU2xvcGUgKyBzcGxpdCA7XG5cbiAgICBzcGxpdCAgPSBEZWdUb1JhZChzcGxpdCApO1xuICAgIHggICA9IHNwcmVhZDBJZHggKihwaGlTbG9wZSArIHBoaVdpbmQqTWF0aC5jb3Moc3BsaXQgKSk7XG4gICAgeSAgID0gc3ByZWFkMElkeCAqKHBoaVdpbmQqTWF0aC5zaW4oc3BsaXQgKSk7XG4gICAgUnYgID0gTWF0aC5zcXJ0KHggKnggICsgeSAqeSApO1xuXG4gICAgc3ByZWFkTWF4ID0gc3ByZWFkMElkeCAgKyBSdiA7XG4gICAgcGhpRXcgPSBzcHJlYWRNYXggLyBzcHJlYWQwSWR4ICAtIDE7XG4gICAgYSAgPSBNYXRoLmFzaW4oTWF0aC5hYnMoeSApIC8gUnYgKTtcbiAgICBpZih4ICA+PSAwLjApXG4gICAgICBhICA9ICh5ICA+PSAwLjApID8gYSAgICAgICAgICAgOiBNX1BJICsgTV9QSSAtIGEgO1xuICAgIGVsc2VcbiAgICAgIGEgID0gKHkgID49IDAuMCkgPyAoTV9QSSAtIGEgKSA6IE1fUEkgKyBhIDtcbiAgICBcbiAgICBzcGxpdCAgPSBSYWRUb0RlZyhhICk7XG4gICAgaWYgKHVwU2xvcGUgKyBzcGxpdCAgPiAzMDYuMClcbiAgICAgIGF6aW11dGhNYXhNYXBbaWR4XSA9IHVwU2xvcGUgKyBzcGxpdCAgLSAzNjAuMDtcbiAgICBlbHNlXG4gICAgICBhemltdXRoTWF4TWFwW2lkeF0gPSB1cFNsb3BlICsgc3BsaXQgO1xuXG4gICAgZWZmZWN0aXZlV2luZCAgPSBNYXRoLnBvdyhwaGlFdypmdWVsUHJvcHMuRnVlbF9XaW5kRSwgMS93aW5kQik7XG4gICAgXG4gICAgLy9EbyBlZmZlY3RpdmUgd2luZCBvbmx5IGlmIHBoaUV3ID4gc21pZGdlblxuICAgIGlmKHBoaUV3ID4gc21pZGdlbikge1xuXG4gICAgICBtYXhXaW5kICA9IDAuOSpyeEludGVuc2l0eU1hcFtpZHhdO1xuICAgICAgaWYoZWZmZWN0aXZlV2luZCAgPiAgbWF4V2luZCApIHtcbiAgICAgICAgcGhpRXcgPSB3aW5kSypNYXRoLnBvdyhtYXhXaW5kICwgd2luZEIpO1xuICAgICAgICBlZmZlY3RpdmVXaW5kICA9IG1heFdpbmQgO1xuICAgICAgICBzcHJlYWRNYXggPSBzcHJlYWQwSWR4ICooMSArIHBoaUV3KTtcbiAgICAgIH1cbiAgICB9XG5cbiAgICBpZihlZmZlY3RpdmVXaW5kICA+ICBzbWlkZ2VuKSB7XG5cbiAgICAgIGx3UmF0aW8gID0gMS4wICsgMC4wMDI4NDA5MDkgKiBlZmZlY3RpdmVXaW5kIDtcbiAgICAgIGlmIChsd1JhdGlvICA+IDEuMDAwMDEpXG4gICAgICAgIGVjY2VudHJpY2l0eU1hcFtpZHhdID0gTWF0aC5zcXJ0KGx3UmF0aW8gKmx3UmF0aW8gIC0gMSkvbHdSYXRpbyA7XG4gICAgfVxuXG4gICAgc3ByZWFkTWF4SWR4ID0gc3ByZWFkTWF4O1xuICAgIHBoaUVmZldpbmRNYXBbaWR4XSA9IHBoaUV3O1xuXG4gIH1cblxuICByZXR1cm4gKCBzcHJlYWRNYXhJZHggKTtcbn1cblxuZnVuY3Rpb24gc3ByZWFkQW55QXppbXV0aChpZHgsIGF6aW11dGgsIHBoaUVmZldpbmRNYXAsIGF6aW11dGhNYXhNYXAsIHJvc01heE1hcCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBlY2NlbnRyaWNpdHlNYXAsIHJvczBNYXAgKVxue1xuXG4gIHZhciBzcHJlYWRBbnk7XG4gIHZhciBlY2NlbnRyaWNpdHk7XG5cbiAgaWYgKHBoaUVmZldpbmRNYXBbaWR4XSA8IHNtaWRnZW4gJiYgYXppbXV0aE1heE1hcFtpZHhdID09PSBhemltdXRoKVxuICBcbiAgICBzcHJlYWRBbnkgPSByb3NNYXhNYXBbaWR4XTtcbiAgXG4gIGVsc2Uge1xuICBcbiAgICBpZiAoKGRpciA9IE1hdGguYWJzKGF6aW11dGhNYXhNYXBbaWR4XSAtIGF6aW11dGgpKSA+IDE4MClcbiAgICAgIGRpciA9IDM2MC4wIC0gZGlyO1xuXG4gICAgZGlyID0gRGVnVG9SYWQoZGlyKTtcblxuICAgIGVjY2VudHJpY2l0eSA9IGVjY2VudHJpY2l0eU1hcFtpZHhdO1xuICAgIHNwcmVhZEFueSA9IHJvc01heE1hcFtpZHhdKigxIC0gZWNjZW50cmljaXR5KS8oMSAtIGVjY2VudHJpY2l0eSpNYXRoLmNvcyhkaXIpKTtcblxuICAgIGlmIChzcHJlYWRBbnkgPiBJTkYpXG4gICAgICBzcHJlYWRBbnkgPSByb3MwTWFwW2lkeF07XG5cbiAgfVxuXG4gIHJldHVybiBzcHJlYWRBbnk7XG59XG5cblxuZnVuY3Rpb24gRGVnVG9SYWQoeCkge1xuICB4ICo9IDAuMDE3NDUzMjkzO1xuICByZXR1cm4geDtcbn1cblxuZnVuY3Rpb24gUmFkVG9EZWcoeCkge1xuICB4ICo9IDU3LjI5NTc3OTUxO1xuICByZXR1cm4geDtcbn1cblxuZnVuY3Rpb24gZXF1YWwoeCx5KXtcbiAgaWYgKCBNYXRoLmFicyh4LXkpPHNtaWRnZW4gKVxuICAgIHJldHVybiB0cnVlO1xuICBlbHNlXG4gICAgcmV0dXJuIGZhbHNlO1xufVxuXG5leHBvcnRzLndpbmRBbmRTbG9wZSA9IHdpbmRBbmRTbG9wZTtcbmV4cG9ydHMuc3ByZWFkQW55QXppbXV0aCA9IHNwcmVhZEFueUF6aW11dGg7XG5leHBvcnRzLm5vV2luZE5vU2xvcGUgPSBub1dpbmROb1Nsb3BlO1xuIl19
;