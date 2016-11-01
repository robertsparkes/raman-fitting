mkdir combinedfigure
mkdir pdf
mkdir jpg
mkdir xy_chart_files

echo "name g_height g_location g_width g_area d1_height d1_location d1_width d1_area d2_height d2_location d2_width d2_area d3_height d3_location d3_width d3_area d4_height d4_location d4_width d4_area r1_ratio r2_ratio r2_temp ra1_ratio ra1_temp ra2_ratio ra2_temp r2voigt plottemp totalwidth totalwidthvoigt" > acombinedresults.txt
echo "name g_height g_location g_full_width g_area d1_height d1_location d1_full-width d1_area d_to_g_intensities" > dgcombinedresults.txt

#makecpt -T300/500/1 > temps.cpt
