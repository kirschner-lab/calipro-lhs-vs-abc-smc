#@ File file_input_roiset
#@ Output tbl

// ROI Manager functions are slow outside of batch mode.
setBatchMode(true);

// Open the image.
open("nejm199501263320402_f1.jpeg");

// Clear any ROIs then loading the file.
if (roiManager("count") > 0) {
    roiManager("deselect");
    roiManager("delete");
}
roiManager("open", file_input_roiset);

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
// Save sa table and CSV.
tbl = File.getNameWithoutExtension(file_input_roiset) + ".csv";
file_output_csv = File.getDirectory(file_input_roiset) + tbl;
Table.create(tbl);
Table.reset(tbl);
Table.setColumn("year", xs);
Table.setColumn("cd4_cells_per_mm3", ys);
Table.update;
Table.save(file_output_csv);

// ROI Manager functions are slow outside of batch mode.
setBatchMode(false);
