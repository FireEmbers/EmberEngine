#Ember Engine: A node module for wildfire simulation

This module is the core function for a  forest fire simulator based on cellular 
automata and fireLib (Rothermel and friends).

##Install
git install git+ssh://git@github.com:FireEmbers/EmberEngine.git

##Usage

To run the module simple do this in a node file:

`var emberEngine = require('emberEngine');`

`emberEngine(dataUnit, Rows, Cols, aspectArray, slopeArray );`

* **Rows** and **Cols** are the size of the mesh, 
* **aspectArray** and **slopeArray** are the terrain aspect and slope
* **dataUnit** is the vector of uniform properties, currently: 
    * dataUnit[0] = FuelMoisture(%)
    * dataUnit[1] = WindVeloctity (m/s) 
    * dataUnit[2] = WindDirection (º clockwise from north) 

Function returns an array of igition times in JSON format

##CrowdProcess usage

The program.js file in src is minified and browserified in build/program.min.js. To create this file just 
run `grunt` in the module directory.

Next, to crete the Run function to pass to the CP API, just do something like this:

```
function RunString(){

  function Run(dataUnit){

    //req is exposed by broserify so that we can have embersEngine module in the Run function
    var engine = req('path_to_embersEngine/src/program.js'); 

    return emberEngine(dataUnit, rows, cols, aspectMap, slopeMap, clcMap, height, width);

  }

  var
   string = Run.toString() + ';' + programString +
  'var rows =' + rows.toString() + ';' +
  'var cols =' + cols.toString() + ';' +
  'var height =' + height.toString() + ';' +
  'var width =' + width.toString() + ';' +
  'var slopeMap =' + JSON.stringify(slopeArray) + ';' +
  'var aspectMap =' + JSON.stringify(aspectArray) + ';' +
  'var clcMap =' + JSON.stringify(clcArray) + ';';

  return string;
}
```

##Test

test tings with 'test' in test folder. eg:

`node test/testFGM.js`

Testing is important. 

In scientific computing testing is called verification. EmberEngine is a port to js of existing
fire models (eg fireLib) so every module change needs to be verified against the original code

