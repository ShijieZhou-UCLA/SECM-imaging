
main.m -- added

clpconfig.m -- added a 7-th parameter "shift" in psf parameters
			   suppressed output

DataSpec.m -- added "start_time", made start_location calculated from start_time
			  added "nskip" as an input parameter that replaces sampleskip
			  added floor(...) in line 65 in definiton of b
			  added zero-stripping of "sampletime" and "currents" in line 89

SecmImage.m -- added method "square_loss(obj, img_truth)" to calculate square loss
			   added method "min_square_loss(obj, x_truth, D)" to calculate minimum square loss up to translation and rotation
			   added method "normalize(obj)" to normalize max intensity to 1
			   added method "translate(obj, xdis, ydis)"

Solver.m -- suppressed display in "display_results(obj)"

IPalm.m -- changed "stop_criterias()" to check every iteration

ReweightIPalm.m -- changed set_niter_ipalm(50) to ..(30), reduced number of ipalm iterations within reweighted ipalm

ReweightIpalmSecmRealdata.m -- supressed display




