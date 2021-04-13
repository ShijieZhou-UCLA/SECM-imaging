**Changelog for Triangle Dictionaries**
---Shijie Zhou 

example6.m -- 

* added. Run this for triangle dictionary real data experiment.

ReweightIPalmMultiMotif.m --

* added. Reweighting for triangle data experiment.

ReweightIPalmSecmRealdataMultiMotif.m --

* added. Reweighting for triangle dictionary real data experiment.  

IPalmSecmSimul.m --
		
* line 31: add "draw_ sum_ images" to display the sum sparse map;
* line 32: add "draw_ images" to display the separate sparse maps;
* line 53 & 54: added only for display separate sparse maps.

ReweightIPalmSecmSimul.m --
		
* line 31: add to display the sum sparse map;
* line 32: add to display the separate sparse maps;
* line 51 & 52: added only for display separate sparse maps.

IPalmSecmRealdata.m --

* line 30: add to display the sum sparse map;
* line 31: add to display the separate sparse maps;
* line 51 & 52: added only for display separate sparse maps.

SecmImageArray.m -- 
						
* line 135: add "isa(obj1, 'cell')" for matrix array for times;
* line 200: add "isa(lda,'cell')" for lambda as a matrix array for soft thresholding;
* line 387: modify function "draw_images";
* line 404: add function "draw_ sum_ images". 

CalibLassoMultiMotif.m --

* line 112: add function "set_ lda_ array".

IPalm.m --

* line 67: add "obj.objval" for smooth part h;

* line 135: add function "set_funchval" to calculate the value of smooth part h.

cell_times.m --

* added. For scalar t times matrix array lda







