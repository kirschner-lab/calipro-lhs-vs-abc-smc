subject_number = 1;

// ROI Manager functions are slow outside of batch mode.
setBatchMode(true);

// Get axis coordinates from the 3 fiducials.
roiManager("select", 0);
Roi.getBounds(x, y, width, height);
origin = newArray(x, y, width, height);

roiManager("select", 1);
Roi.getBounds(x, y, width, height);
x_end = newArray(x, y, width, height);

roiManager("select", 2);
Roi.getBounds(x, y, width, height);
y_end = newArray(x, y, width, height);

// Axis calibrations.
x_len = x_end[0] - origin[0];
x_cal = 10 / x_len;

y_len = origin[1] - y_end[1];
y_cal = 2000 /y_len;

// Number of points (minus 3 for axis fiducials).
n = roiManager("count") - 3;
xs = newArray(n);
ys = newArray(n);
for(i = 0; i < n; i++) {
	roiManager("select", i + 3);
	Roi.getBounds(x, y, width, height);
	xs[i] = ((x + width/2) - (origin[0] + origin[2]/2)) * x_cal;
	ys[i] = ((origin[1] + origin[3]/2) - (y + height/2)) * y_cal;
}
tbl = "subject-" + subject_number;
Table.create(tbl);
Table.reset(tbl);
Table.setColumn("x", xs);
Table.setColumn("y", ys);
Table.update;

// ROI Manager functions are slow outside of batch mode.
setBatchMode(false);