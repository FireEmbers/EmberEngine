
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
var fireLib = require('./fireLib');

module.exports = function (dataArray, rowsPC, colsPC, aspectMapPC, slopeMapPC, clcMapPC, heightPC, widthPC){

  ignIdx =  Math.floor(colsPC/2) + Math.floor(rowsPC/2)*colsPC;

  if (!(/\b32\d\b/.test(clcMapPC[ignIdx]) || /\b31\d\b/.test(clcMapPC[ignIdx])) ){
      return null;
  }

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
    ignMap[ignIdx] = 0;



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