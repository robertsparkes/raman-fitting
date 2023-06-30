#!/bin/bash

###################################### Description #################################################
#
# This script analyses Raman spectra of Carbonaceous Material, by fitting Lorentzian distributions
# to the G, D1, D2, D3 and D4 peaks, as well as correcting for a linear background.
# The input is taken from a series of text files as produced by a Renishaw Raman spectrometer using
# Wire software. Proprietary .wxd files can be converted into two-column space-separated text files 
# (wavenumber intensity) using the "Wire Batch Convert" program. The text files should be contained 
# within one single folder, or grouped into sub-folders.
#
# The script outputs three graphs, containing the raw spectra with linear background identified,
# raw spectra with overalll fit superimposed and a residual shown, and the spectra following the fitting,
# showing the fitted peaks after the background has been removed. The fitting parameters (peak locations,
# amplitudes, widths and areas, as well as characteristic area ratios) are outputted to a summary file 
# for further analysis
#
# The script requires the following software to run:
# - A Unix / Linux environment (tested with Ubuntu)
# - Bash terminal program
# - Dos2unix text file conversion software
# - Gnuplot graphing software. 
#     - Version 4.5 or above is required
# - Ghostscript PostScript and PDF manipulation software
# - The script "prepraman.sh" should be run before the first files in a given folder are analysed.
#   This script creates folders and initiates some datafiles for the subsequent fits
# - Both this fitting script and "prepraman.sh" require permission to execute as programs
#
# The script executes from the command line, in the form
# $ sparkesfitraman.sh [options] [input files]
#
# The options are
# -q Quiet mode - graphs appear on screen but immediately disappear
# -d Delete - removes previous files from "acombinedresults.txt"
# -t[value] Threshold - the signal-to-noise ratio below which a peak is too noisy to process
#
# Input files can be listed individually, or selected all at once using a wildcard (e.g. *.txt)
# After analysis the results are written to a file entitled "acombinedresults.txt". Any filename
# already in this file will be ignored and not re-fitted, hence the delete option.
#
# Example code to prepare for and then analyse all samples with "taiwan" in the file name:
# $ prepraman.sh 
# $ sparkesfitraman.sh -d -q -t 5 taiwan*.txt
#
# Note: The included script "cropraman.sh" will take a file and crop to certain wavenumbers. The
# fitting procedure is less accurate if files extend too far beyond 1900 cm-1 as the assumption
# of a linear background is no longer valid.
#
##################################################################################################


# First the options are collected


quiet=false
persist=-persist
delete=false
delans=n
threshold=2
silthreshold=2

# r2 limit should be 100 times greater (and an integer)
r2limit=60

while getopts 'dqt:' option
do
case $option in
	d) delete=true;;
	q) quiet=true;;
	t) threshold=$OPTARG;;
esac
done
shift $(($OPTIND - 1))

if [ "$delete" = "true" ] ; then
	echo "Really delete all records? (y/n)"
	read delans
	if [ "$delans" = "y" ] ; then
		echo "name g_height g_location g_width g_area d1_height d1_location d1_width d1_area d2_height d2_location d2_width d2_area d3_height d3_location d3_width d3_area d4_height d4_location d4_width d4_area r1_ratio r2_ratio r2_temp ra1_ratio ra1_temp ra2_ratio ra2_temp r2voigt plottemp totalwidth totalwidthvoigt" > acombinedresults.txt

		echo "Records deleted!"
	else
		echo "Records saved!"	
	fi
fi

echo "This file reports in FWHM" >> acombinedresults.txt

echo Quiet? $quiet

if [ "$quiet" = "true" ] ; then
persist=""
fi

echo $#
echo $@

function processsample {
#Prepare input file
filename=$1
dos2unix $filename
nicename=${filename%\.*}
echo $nicename
shortname=`echo $nicename | awk 'BEGIN { FS = "_" } { print $2 }'`
echo $shortname

outputresult=`awk 'match($1, '/$nicename/')' acombinedresults.txt`
if [ "$outputresult" = "" ]; 
then 

rm fit.log






########################    Read basic spectrum parameters     ##########################
yinit2=`awk 'END {print $2}' $1`
yend=`awk 'NR==1 {print $2}' $1`
xinit=`awk 'END {print $1}' $1`
xend=`awk 'NR==1 {print $1}' $1`

grad=`echo "scale=2; ($yend - $yinit2)/($xend - $xinit)" | bc`

yinit=`echo "scale=2; $yinit2 - ($grad * $xinit)" | bc`

echo background = $grad x + $yinit




###############      Calculate Signal to Noise Ratio, exit if too noisy  ################

snrmax=`awk 'BEGIN { max = 1 } (($1 > 1740) && ($1 < 1830)) { if ( max < $2 ) max = $2 } END { print max }' $1`
snrmin=`awk 'BEGIN { min = '$snrmax' } (($1 > 1740) && ($1 < 1830)) { if ( min > $2 ) min = $2 } END { print min - 0.1 }' $1`
signalmaxy=`awk 'BEGIN { max = 0 } (( $1 > 1200 ) && ( $1 < 1790 )) { if ( max < $2 ) max = $2 } END { print max }' $1`
signalmaxx=`awk ' $2 ~ /'$signalmaxy'/ { print $1 }' $1`
signalmax=`awk 'BEGIN { max = 0 } (( $1 > 1200 ) && ( $1 < 1790 )) { if ( max < $2 ) max = $2 } END { print max - ( '$signalmaxx' * '$grad' + '$yinit' ) }' $1`

snr=`echo "scale=0; $signalmax / ($snrmax - $snrmin)" | bc`
echo Noise highpoint snrmax = $snrmax
echo Noise lowpoint snrmin = $snrmin
echo Signal peak location signalmaxx = $signalmaxx
echo Uncorrected signal peak signalmaxy = $signalmaxy
echo Corrected signal peak signalmax = $signalmax
echo snr = $snr

if (( "$snr" < "$threshold" ))
then

echo "Noisy / no signal"

gnuplot <<EOF

set term post landscape color solid 8
set output 'combined.ps'

set title "$nicename = Noisy"

plot '$1' with lines title "$nicename = Noisy"



EOF

ps2pdf combined.ps ${nicename}combined.pdf
rm combined.ps

echo $nicename Noisy >> acombinedresults.txt

return

fi


############################ GNUPlot Curve Fitting  ###########################################

d1height=`awk 'BEGIN { max = -100000000000 } (($1 > 1200) && ($1 < 1450)) { if ( max < $2 ) max = $2 } END { print 10 * ( max - ( '$yinit' + 1350 * '$grad' ) ) }' $1`
d1heightb=`awk 'BEGIN { max = -100000000000 } (($1 > 1200) && ($1 < 1450)) { if ( max < $2 ) max = $2 } END { print 15 * max }' $1`
d1loc=`awk ' /'$d1heightb'/ && ( $1 > 1200 ) && ( $1 < 1450 ) { print $1 }' $1`
d2height=`awk 'BEGIN { max = -100000000000 } (($1 > 1605) && ($1 < 1640)) { if ( max < $2 ) max = $2 } END { print 10 * ( max - ( '$yinit' + 1600 * '$grad' ) ) }' $1`
d3height=`awk 'BEGIN { max = -100000000000 } (($1 > 1490) && ($1 < 1510)) { if ( max < $2 ) max = $2 } END { print 10 * ( max - ( '$yinit'+ 1500 * '$grad' ) ) }' $1`
gheight=`awk 'BEGIN { max = -100000000000 } (($1 > 1575) && ($1 < 1600)) { if ( max < $2 ) max = $2 } END { print max - ( '$yinit' + 1600 * '$grad' ) }' $1`
gheightb=`awk 'BEGIN { max = -100000000000 } (($1 > 1575) && ($1 < 1600)) { if ( max < $2 ) max = $2 } END { print 10 * max }' $1`
gloc=`awk ' /'$gheightb'/  && ( $1 > 1550 ) && ( $1 < 1650 ){ print $1 }' $1`
#d2height=`awk 'BEGIN { max = -100000000000 } (($1 > 1610) && ($1 < 1620)) { if ( max < $2 ) max = $2 } END { print max - '$yinit' - '$gheight' }' $1`

# Remember to change these lines if d1 restrictions change
sind1loc=`echo "s ( ( ( 3.14159 * ( $d1loc - 1300 ) / 100 ) - ( 3.14159 / 2 ) ) )" | bc -l`
cosd1loc=`echo "c ( ( ( 3.14159 * ( $d1loc - 1300 ) / 100 ) - ( 3.14159 / 2 ) ) )" | bc -l`
tand1loc=`echo "scale=5; $sind1loc / $cosd1loc" | bc`

# Remember to change these lines if g restrictions change
singloc=`echo "s ( ( ( 3.14159 * ( $gloc - 1563 ) / 60 ) - ( 3.14159 / 2 ) ) )" | bc -l`
cosgloc=`echo "c ( ( ( 3.14159 * ( $gloc - 1563 ) / 60 ) - ( 3.14159 / 2 ) ) )" | bc -l`
tangloc=`echo "scale=5; $singloc / $cosgloc" | bc`


echo d1height = $d1height
echo d1loc = $d1loc = $tand1loc
echo d2height = $d2height
echo gheight = $gheight
echo gloc = $gloc

echo grad = $grad > param.txt
echo int = $yinit >> param.txt


echo gheight = $gheightb >> param.txt
echo gwidth = -5 >> param.txt
echo gloc = $tangloc >> param.txt

echo d1height = $d1height >> param.txt
echo d1width = -5 >> param.txt
echo d1loc = 0.1 >> param.txt


echo d2height = $d2height >> param.txt
echo d2width = -5 >> param.txt
echo d2loc = 0.6 >> param.txt

echo d3loc = 0.1 >> param.txt
echo d3height = $d3height >> param.txt
echo d3width = -5 >> param.txt




gnuplot $persist<<EOF

# Restrict gloc to the range of [1563:1605] and amplitude positive with width up to 80 cm-1
g(x) = (1605-1563)/pi*(atan(x)+pi/2)+1563
gh(x) = sqrt(x**2)
gw(x) = 40/pi*(atan(x)+pi/2)+0.1

pr g($tangloc)

# Restrict d2 to the range of [1605:1625] and width to range 1-32cm
d2(x) = (1625-1605)/pi*(atan(x)+pi/2)+1605
d2h(x) = sqrt(x**2)
d2w(x) = 16/pi*(atan(x)+pi/2)+0.1

# Restrict d1 to the range of [1300:1400] and width to range 1-200cm
d1(x) = (1400-1300)/pi*(atan(x)+pi/2)+1300
d1h(x) = sqrt(x**2)
d1w(x) = 100/pi*(atan(x)+pi/2)+0.1

# Restrict d3 to the range of [1475-1525] and width to range 1-200cm
d3(x) = (1525-1475)/pi*(atan(x)+pi/2)+1475
d3h(x) = sqrt(x**2)
d3w(x) = 100/pi*(atan(x)+pi/2)+1

bg(x) = int + grad * x

gpeak(x) = gh(gheight) * voigt( x - g(gloc) , gw(gwidth) )

d1peak(x) = d1h(d1height) * voigt( x - d1(d1loc) , d1w(d1width) )

d2peak(x) = d2h(d2height) * voigt( x - d2(d2loc) , d2w(d2width) )

#d3peak(x) = d3h(d3height) * voigt( x - d3(d3loc) , d3w(d3width) )

d3peak(x) = 0 * x
d4peak(x) = 0 * x

# f(x) = gpeak(x) + d1peak(x) + d2peak(x) + d3peak(x) + bg(x)
f(x) = gpeak(x) + d1peak(x) + d2peak(x) + bg(x)

p(x) = gpeak(x) + d1peak(x) + d2peak(x)

#plot [x=1000:1900] '$1' using 1:(\$2 - bg(\$1)), p(x)

FIT_LIMIT = 1e-9
FIT_MAXITER = 1000
fit f(x) '$1' using 1:2:(1) via 'param.txt'

set table "residual.xy"
plot [x=1000:1900] '$1' using 1:(\$2 - f(\$1))

set table "lorentzians.xy"
plot [x=1000:1900] p(x)

set table "bgremoved.xy"
plot [x=1000:1900] '$1' using 1:(\$2 - bg(\$1))

set table "d1peak.xy"
plot [x=1000:1900] d1peak(x)

set table "d2peak.xy"
plot [x=1000:1900] d2peak(x)

set table "gpeak.xy"
plot [x=1000:1900] gpeak(x)

set table "bgremoved.xy"
plot [x=1000:1900] '$1' using 1:(\$2 - bg(\$1))

pr "gloc ", g(gloc)
pr "garea ", gh(gheight)
pr "gwidth ", gw(gwidth)
pr "gheight ", gh(gheight)*voigt(0,gw(gwidth))

pr "d1loc ", d1(d1loc)
pr "d1area ", d1h(d1height)
pr "d1width ", d1w(d1width)
pr "d1height ", d1h(d1height)*voigt(0,d1w(d1width))

pr "d2loc ", d2(d2loc)
pr "d2area ", d2h(d2height)
pr "d2width ", d2w(d2width)
pr "d2height ", d2h(d2height)*voigt(0,d2w(d2width))

pr "d3loc ", d3(d3loc)
pr "d3area ", d3h(d3height)
pr "d3width ", d3w(d3width)
pr "d3height ", d3h(d3height)*voigt(0,d3w(d3width))

pr "int ", int
pr "grad ", grad


set print 'param3a.txt'
pr "gloc ", g(gloc)
pr "garea ", gh(gheight)
pr "gwidth ", gw(gwidth)
pr "gheight ", gh(gheight)*voigt(0,gw(gwidth))

pr "d1loc ", d1(d1loc)
pr "d1area ", d1h(d1height)
pr "d1width ", d1w(d1width)
pr "d1height ", d1h(d1height)*voigt(0,d1w(d1width))

pr "d2loc ", d2(d2loc)
pr "d2area ", d2h(d2height)
pr "d2width ", d2w(d2width)
pr "d2height ", d2h(d2height)*voigt(0,d2w(d2width))

pr "d3loc ", d3(d3loc)
pr "d3area ", d3h(d3height)
pr "d3width ", d3w(d3width)
pr "d3height ", d3h(d3height)*voigt(0,d3w(d3width))

pr "int ", int
pr "grad ", grad


save fit 'param.txt' 'param2.txt'

save "savefilevoigt.plt"

EOF

#Remove scientific notation
sed 's/e-/\*10\^-/' param3a.txt > param3.txt

# Find final parameters #
gloc=`awk ' $1 ~ /gloc/ { print sqrt( $2 ^ 2 ) } ' param3.txt `
gheight=`awk ' $1 ~ /gheight/ { print $2 } ' param3.txt `
gwidth=`awk ' ( NR > 4) $1 ~ /gwidth/ { print $2 } ' param3.txt `
garea=`awk ' ( NR > 4) $1 ~ /garea/ { print $2 } ' param3.txt `


d1loc=`awk ' $1 ~ /d1loc/ { print sqrt( $2 ^ 2 ) } ' param3.txt `
d1height=`awk ' $1 ~ /d1height/ { print $2 } ' param3.txt `
d1width=`awk ' $1 ~ /d1width/ { print $2 } ' param3.txt `
d1area=`awk ' ( NR > 4) $1 ~ /d1area/ { print $2 } ' param3.txt `

d2loc=`awk ' $1 ~ /d2loc/ { print sqrt( $2 ^ 2 ) } ' param3.txt `
d2height=`awk ' $1 ~ /d2height/ { print $2 } ' param3.txt `
d2width=`awk ' $1 ~ /d2width/ { print $2 } ' param3.txt `
d2area=`awk ' ( NR > 4) $1 ~ /d2area/ { print $2 } ' param3.txt `

d3loc=`awk ' $1 ~ /d3loc/ { print sqrt( $2 ^ 2 ) } ' param3.txt `
d3height=`awk ' $1 ~ /d3height/ { print $2 } ' param3.txt `
d3width=`awk ' $1 ~ /d3width/ { print $2 } ' param3.txt `
d3area=`awk ' ( NR > 4) $1 ~ /d3area/ { print $2 } ' param3.txt `

bggrad=`awk ' $1 ~ /grad/ { print $2 } ' param3.txt `
bgint=`awk ' $1 ~ /int/ { print $2 } ' param3.txt `

echo G height = $gheight
echo G location = $gloc
echo G width = $gwidth
echo G area = $garea

echo D1 height == area = $d1height
echo D1 location = $d1loc
echo D1 width = $d1width
echo D1 area = $d1area

echo D2 height == area = $d2height
echo D2 location = $d2loc
echo D2 width = $d2width
echo D2 area = $d2area

echo Lorentz solid offset = $bgint
echo Lorentz inclined offset = $bggrad

#Calculate SNR

snrmax=`awk 'BEGIN { max = -100000000000 } (($1 > 1700) && ($1 < 1800)) { if ( max < $2 ) max = $2 } END { print max }' bgremoved.xy`
snrmin=`awk 'BEGIN { min = 100000000000 } (($1 > 1700) && ($1 < 1800)) { if ( min > $2 ) min = $2 } END { print min }' bgremoved.xy`
bgremovedmax=`awk 'BEGIN { max = -100000000000 } ( $1 > 1200 ) { if ( max < $2 ) max = $2 } END { print max }' bgremoved.xy`

snr=`echo "scale=0; $bgremovedmax / ($snrmax - $snrmin)" | bc`
echo snrmax = $snrmax
echo snrmin = $snrmin
echo bgremovedmax = $bgremovedmax
echo snr = $snr

# R1 ratio
r1ratio=`echo "scale=5; $d1height / $gheight" | bc`
r1ratiocheck=`echo "scale=0; 100 * $r1ratio / 1" | bc`
echo R1 ratio = $r1ratio $r1ratiocheck

r2ratio=`echo "scale=5; $d1area / ( $d1area + $garea + $d2area )" | bc`
r2ratioa=`echo "scale=2; $d1area / ( $d1area + $garea + $d2area )" | bc`
r2ratiovoigt=$r2ratio
totalwidthvoigt=`echo "scale=2; $gwidth + $d1width + $d2width" | bc`

r2ratiocheck=`echo "scale=0; 100 * $r2ratioa / 1" | bc`
echo R2 ratio = $r2ratio
r2tempa=`echo "scale=3; ((-445 * $r2ratio) + 641 ) / 1 " | bc`
echo R2 temp = $r2tempa

d1widthcheck=`echo "scale=0; $d1width / 1" | bc`


if (( "$r2ratiocheck" < "$r2limit" && "$d1widthcheck" < "60" ))
then
fitstyle=Voigt
temp="Standard R2 temp $r2tempa"
plottemp=$r2tempa

ra1ratio=na
ra2ratio=na
ra1temp=na
ra2temp=na
d3loc=na
d3height=na
d3width=na
d3area=na
d4loc=na
d4height=na
d4width=na
d4area=na

mv param3a.txt paramvoigt.txt

totalwidth=`echo "scale=2; $gwidth + $d1width + $d2width" | bc`

savefile=savefilevoigt.plt

elif (( "$r2ratiocheck" < "$r2limit" && "$r1ratiocheck" < "50" ))
then

fitstyle=Voigt3
temp="Standard R2 temp $r2tempa"
plottemp=$r2tempa

ra1ratio=na
ra2ratio=na
ra1temp=na
ra2temp=na
d3loc=na
d3height=na
d3width=na
d3area=na
d4loc=na
d4height=na
d4width=na
d4area=na

mv param3a.txt paramvoigt.txt

totalwidth=`echo "scale=2; $gwidth + $d1width + $d2width" | bc`

savefile=savefilevoigt.plt

else

echo "R2 or D1width too high for Voigt fit, trying Lorentzians"

mv param3a.txt paramvoigt.txt

######################### Fit using Lorentzians  ###############################
rm fit.log

########################    GNUPLOT Curve Fitting     ##########################
yinit2=`awk 'END {print $2}' $1`
yend=`awk 'NR==1 {print $2}' $1`
xinit=`awk 'END {print $1}' $1`
xend=`awk 'NR==1 {print $1}' $1`

grad=`echo "scale=2; ($yend - $yinit2)/($xend - $xinit)" | bc`

yinit=`echo "scale=2; $yinit2 - ($grad * $xinit)" | bc`

echo background = $grad x + $yinit

d1height=`awk 'BEGIN { max = -100000000000 } (($1 > 1200) && ($1 < 1450)) { if ( max < $2 ) max = $2 } END { print max - ( '$yinit' + 1350 * '$grad' ) }' $1`
d1heightb=`awk 'BEGIN { max = -100000000000 } (($1 > 1200) && ($1 < 1450)) { if ( max < $2 ) max = $2 } END { print max }' $1`
d1loc=`awk ' /'$d1heightb'/ && ( $1 > 1200 ) && ( $1 < 1450 ) { print $1 }' $1`
d2height=`awk 'BEGIN { max = -100000000000 } (($1 > 1610) && ($1 < 1640)) { if ( max < $2 ) max = $2 } END { print max - ( '$yinit' + 1600 * '$grad' ) }' $1`
d3height=`awk 'BEGIN { max = -100000000000 } (($1 > 1490) && ($1 < 1510)) { if ( max < $2 ) max = $2 } END { print max - ( '$yinit'+ 1500 * '$grad' ) }' $1`
d4height=`awk 'BEGIN { max = -100000000000 } (($1 > 1140) && ($1 < 1150)) { if ( max < $2 ) max = $2 } END { print max - ( '$yinit' + 1150 * '$grad' ) }' $1`
gheight=`awk 'BEGIN { max = -100000000000 } (($1 > 1580) && ($1 < 1600)) { if ( max < $2 ) max = $2 } END { print max - ( '$yinit' + 1600 * '$grad' ) }' $1`
gheightb=`awk 'BEGIN { max = -100000000000 } (($1 > 1580) && ($1 < 1600)) { if ( max < $2 ) max = $2 } END { print max }' $1`
gloc=`awk ' /'$gheightb'/  && ( $1 > 1550 ) && ( $1 < 1650 ){ print $1 }' $1`
#d2height=`awk 'BEGIN { max = -100000000000 } (($1 > 1610) && ($1 < 1620)) { if ( max < $2 ) max = $2 } END { print max - '$yinit' - '$gheight' }' $1`

# Remember to change these lines if d1 restrictions change
sind1loc=`echo "s ( ( ( 3.14159 * ( $d1loc - 1300 ) / 100 ) - ( 3.14159 / 2 ) ) )" | bc -l`
cosd1loc=`echo "c ( ( ( 3.14159 * ( $d1loc - 1300 ) / 100 ) - ( 3.14159 / 2 ) ) )" | bc -l`
tand1loc=`echo "scale=5; $sind1loc / $cosd1loc" | bc`

# Remember to change these lines if g restrictions change
singloc=`echo "s ( ( ( 3.14159 * ( $gloc - 1567 ) / 40 ) - ( 3.14159 / 2 ) ) )" | bc -l`
cosgloc=`echo "c ( ( ( 3.14159 * ( $gloc - 1567 ) / 40 ) - ( 3.14159 / 2 ) ) )" | bc -l`
tangloc=`echo "scale=5; $singloc / $cosgloc" | bc`


echo d1height = $d1height
echo d1loc = $d1loc = $tand1loc
echo d2height = $d2height
echo gheight = $gheight
echo gloc = $gloc

echo grad = $grad > param.txt
echo int = $yinit >> param.txt


echo gheight = $gheight >> param.txt
echo gwidth = -1.5 >> param.txt
echo gloc = $tangloc >> param.txt

echo d1height = $d1height >> param.txt
echo d1width = -0.5 >> param.txt
echo d1loc = $tand1loc >> param.txt


echo d2height = $d2height >> param.txt
echo d2width = -1.5 >> param.txt
echo d2loc = -5 >> param.txt

echo d3loc = 0.1 >> param.txt
echo d3height = $d3height >> param.txt
echo d3width = 1 >> param.txt

echo d4loc = 5 >> param.txt
echo d4height = $d4height >> param.txt
echo d4width = 1 >> param.txt





gnuplot $persist<<EOF

# Restrict gloc to the range of [1567:1605] and amplitude positive with width up to 80 cm-1
g(x) = (1605-1567)/pi*(atan(x)+pi/2)+1567
gh(x) = sqrt(x**2)
gw(x) = 40/pi*(atan(x)+pi/2)+1

pr g($tangloc)

# Restrict d2 to the range of [1600:1630] and width to range 1-80cm
d2(x) = (1630-1600)/pi*(atan(x)+pi/2)+1600
d2h(x) = sqrt(x**2)
d2w(x) = 40/pi*(atan(x)+pi/2)+1

# Restrict d1 to the range of [1300:1400] and width to range 1-200cm
d1(x) = (1400-1300)/pi*(atan(x)+pi/2)+1300
d1h(x) = sqrt(x**2)
d1w(x) = 100/pi*(atan(x)+pi/2)+1

# Restrict d3 to the range of [1475-1525] and width to range 1-200cm
d3(x) = (1525-1475)/pi*(atan(x)+pi/2)+1475
d3h(x) = sqrt(x**2)
d3w(x) = 100/pi*(atan(x)+pi/2)+1

# Restrict d4 to the range of [1200-1250] and width to range 1-200cm
d4(x) = (1250-1200)/pi*(atan(x)+pi/2)+1200
d4h(x) = sqrt(x**2)
d4w(x) = 100/pi*(atan(x)+pi/2)+1

bg(x) = int + grad * x

gpeak(x) = gh(gheight) * ( (gw(gwidth))**2 / (( x - g(gloc))**2 + (gw(gwidth))**2 ))

d1peak(x) = d1h(d1height) * ( (d1w(d1width))**2 / (( x - d1(d1loc))**2 + (d1w(d1width))**2 ))

d2peak(x) = d2h(d2height) * ( (d2w(d2width))**2 / (( x - d2(d2loc))**2 + (d2w(d2width))**2 ))

d3peak(x) = d3h(d3height) * ( d3w((d3width))**2 / (( x - d3(d3loc) )**2 + (d3w(d3width))**2 ))

d4peak(x) = d4h(d4height) * ( d4w((d4width))**2 / (( x - d4(d4loc) )**2 + (d4w(d4width))**2 ))

f(x) = gpeak(x) + d1peak(x) + d2peak(x) + d3peak(x) + d4peak(x) + bg(x)

p(x) = gpeak(x) + d1peak(x) + d2peak(x) + d3peak(x) + d4peak(x)

FIT_LIMIT = 1e-9
FIT_MAXITER = 1000
fit f(x) '$1' using 1:2:(1) via 'param.txt'

set table "residual.xy"
plot [x=800:2200] '$1' using 1:(\$2 - f(\$1))

set table "lorentzians.xy"
plot [x=800:2200] p(x)

set table "bgremoved.xy"
plot [x=800:2200] '$1' using 1:(\$2 - bg(\$1))

set table "d1peak.xy"
plot [x=800:2200] d1peak(x)

set table "d2peak.xy"
plot [x=800:2200] d2peak(x)

set table "gpeak.xy"
plot [x=800:2200] gpeak(x)

set table "bgremoved.xy"
plot [x=800:2200] '$1' using 1:(\$2 - bg(\$1))

pr "gloc ", g(gloc)
pr "gheight ", gh(gheight)
pr "gwidth ", gw(gwidth)

pr "d1loc ", d1(d1loc)
pr "d1height", d1h(d1height)
pr "d1width ", d1w(d1width)

pr "d2loc ", d2(d2loc)
pr "d2height", d2h(d2height)
pr "d2width ", d2w(d2width)

pr "d3loc ", d3(d3loc)
pr "d3height", d3h(d3height)
pr "d3width ", d3w(d3width)

pr "d4loc ", d4(d4loc)
pr "d4height", d4h(d4height)
pr "d4width ", d4w(d4width)

pr "int", int
pr "grad", grad


set print 'param3.txt'
pr "gloc ", g(gloc)
pr "gheight ", gh(gheight)
pr "gwidth ", gw(gwidth)

pr "d1loc ", d1(d1loc)
pr "d1height ", d1h(d1height)
pr "d1width ", d1w(d1width)

pr "d2loc ", d2(d2loc)
pr "d2height ", d2h(d2height)
pr "d2width ", d2w(d2width)

pr "d3loc ", d3(d3loc)
pr "d3height ", d3h(d3height)
pr "d3width ", d3w(d3width)

pr "d4loc ", d4(d4loc)
pr "d4height ", d4h(d4height)
pr "d4width ", d4w(d4width)

pr "int ", int
pr "grad ", grad


save fit 'param.txt' 'param2.txt'

save "savefilelor.plt"

EOF

#Remove scientific notation
#sed 's/e-/\*10\^-/' param3a.txt > param3.txt

# Find final parameters #
gloc=`awk ' $1 ~ /gloc/ { print $2 } ' param3.txt `
gheight=`awk ' $1 ~ /gheight/ { printf "%i", $2 } ' param3.txt `
gwidth=`awk ' $1 ~ /gwidth/ { print $2 } ' param3.txt `

d1loc=`awk ' $1 ~ /d1loc/ { printf "%i", sqrt( $2 ^ 2 ) } ' param3.txt `
d1height=`awk ' $1 ~ /d1height/ { printf "%i", $2 } ' param3.txt `
d1width=`awk ' $1 ~ /d1width/ { print $2 } ' param3.txt `

d2loc=`awk ' $1 ~ /d2loc/ { printf "%i", sqrt( $2 ^ 2 ) } ' param3.txt `
d2height=`awk ' $1 ~ /d2height/ { printf "%i", $2 } ' param3.txt `
d2width=`awk ' $1 ~ /d2width/ { print $2 } ' param3.txt `

d3loc=`awk ' $1 ~ /d3loc/ { printf "%i", sqrt( $2 ^ 2 ) } ' param3.txt `
d3height=`awk ' $1 ~ /d3height/ { printf "%i", $2 } ' param3.txt `
d3width=`awk ' $1 ~ /d3width/ { print $2 } ' param3.txt `

d4loc=`awk ' $1 ~ /d4loc/ { printf "%i", sqrt( $2 ^ 2 ) } ' param3.txt `
d4height=`awk ' $1 ~ /d4height/ { printf "%i", $2 } ' param3.txt `
d4width=`awk ' $1 ~ /d4width/ { print $2 } ' param3.txt `

bggrad=`awk ' $1 ~ /grad/ { print $2 } ' param3.txt `
bgint=`awk ' $1 ~ /int/ { print $2 } ' param3.txt `

echo Lorentz g mean = $gloc
echo Lorentz g amplitude = $gheight
echo Lorentz g width = $gwidth

echo Lorentz d1 mean = $d1loc
echo Lorentz d1 amplitude = $d1height
echo Lorentz d1 width = $d1width

echo Lorentz d2 mean = $d2loc
echo Lorentz d2 amplitude = $d2height
echo Lorentz d2 width = $d2width

echo Lorentz solid offset = $bgint
echo Lorentz inclined offset = $bggrad

#Calculate SNR

snrmax=`awk 'BEGIN { max = -100000000000 } (($1 > 1700) && ($1 < 1800)) { if ( max < $2 ) max = $2 } END { print max }' bgremoved.xy`
snrmin=`awk 'BEGIN { min = 100000000000 } (($1 > 1700) && ($1 < 1800)) { if ( min > $2 ) min = $2 } END { print min }' bgremoved.xy`
bgremovedmax=`awk 'BEGIN { max = -100000000000 } ( $1 > 1200 ) { if ( max < $2 ) max = $2 } END { print max }' bgremoved.xy`

snr=`echo "scale=0; $bgremovedmax / ($snrmax - $snrmin)" | bc`
echo snrmax = $snrmax
echo snrmin = $snrmin
echo bgremovedmax = $bgremovedmax
echo snr = $snr


#Calculate Areas
garea=`echo "scale=5; $gheight * 3.14159 * $gwidth" | bc`
d1area=`echo "scale=5; $d1height * 3.14159 * $d1width" | bc`
d2area=`echo "scale=5; $d2height * 3.14159 * $d2width" | bc`
d3area=`echo "scale=5; $d3height * 3.14159 * $d3width" | bc`
d4area=`echo "scale=5; $d4height * 3.14159 * $d4width" | bc`

echo G area = $garea
echo D1 area = $d1area
echo D2 area = $d2area
echo D3 area = $d3area
echo D4 area = $d4area

r1ratio=`echo "scale=5; $d1height / $gheight" | bc`

r2ratio=`echo "scale=5; $d1area / ( $d1area + $garea + $d2area )" | bc`
r2ratioa=`echo "scale=2; $d1area / ( $d1area + $garea + $d2area )" | bc`
r2ratiob=`echo "scale=0; $r2ratioa * 1000 / 1 " | bc`
echo R2 ratio = $r2ratio
r2tempa=`echo "scale=3; ((-445 * $r2ratio) + 641 ) / 1 " | bc`
echo R2 temp = $r2tempa

ra1ratio=`echo "scale=5; ( $d1area + $d4area ) / ( $d1area + $garea + $d2area + $d3area +$d4area )" | bc`
ra1ratioa=`echo "scale=3; ( $d1area + $d4area ) / ( $d1area + $garea + $d2area + $d3area +$d4area )" | bc`
echo RA1 ratio = $ra1ratioa
ra1temp=`echo "scale=3; ( $ra1ratio - 0.3758 ) / 0.0008 " | bc`
echo RA1 temp = $ra1temp

ra2ratio=`echo "scale=5; ( $d1area + $d4area ) / ( $garea + $d2area + $d3area )" | bc`
ra2ratioa=`echo "scale=3; ( $d1area + $d4area ) / ( $garea + $d2area + $d3area )" | bc`
echo RA2 ratio = $ra2ratioa
ra2temp=`echo "scale=3; ( $ra2ratio - 0.27 ) / 0.0045 " | bc`
echo RA2 temp = $ra2temp
ra2ratiocheck=`echo "scale=0; $ra2ratio * 100 / 1 " | bc`


if (( "$ra2ratiocheck" > "200" ))
then

echo "RA2 too high, returning to Voigt"

# Find Voigt parameters #
gloc=`awk ' $1 ~ /gloc/ { print sqrt( $2 ^ 2 ) } ' paramvoigt.txt `
gheight=`awk ' $1 ~ /gheight/ { print $2 } ' paramvoigt.txt `
gwidth=`awk ' ( NR > 4) $1 ~ /gwidth/ { print $2 } ' paramvoigt.txt `
garea=`awk ' ( NR > 4) $1 ~ /garea/ { print $2 } ' paramvoigt.txt `


d1loc=`awk ' $1 ~ /d1loc/ { print sqrt( $2 ^ 2 ) } ' paramvoigt.txt `
d1height=`awk ' $1 ~ /d1height/ { print $2 } ' paramvoigt.txt `
d1width=`awk ' $1 ~ /d1width/ { print $2 } ' paramvoigt.txt `
d1area=`awk ' ( NR > 4) $1 ~ /d1area/ { print $2 } ' paramvoigt.txt `

d2loc=`awk ' $1 ~ /d2loc/ { print sqrt( $2 ^ 2 ) } ' paramvoigt.txt `
d2height=`awk ' $1 ~ /d2height/ { print $2 } ' paramvoigt.txt `
d2width=`awk ' $1 ~ /d2width/ { print $2 } ' paramvoigt.txt `
d2area=`awk ' ( NR > 4) $1 ~ /d2area/ { print $2 } ' paramvoigt.txt `

d3loc=`awk ' $1 ~ /d3loc/ { print sqrt( $2 ^ 2 ) } ' paramvoigt.txt `
d3height=`awk ' $1 ~ /d3height/ { print $2 } ' paramvoigt.txt `
d3width=`awk ' $1 ~ /d3width/ { print $2 } ' paramvoigt.txt `
d3area=`awk ' ( NR > 4) $1 ~ /d3area/ { print $2 } ' paramvoigt.txt `

bggrad=`awk ' $1 ~ /grad/ { print $2 } ' paramvoigt.txt `
bgint=`awk ' $1 ~ /int/ { print $2 } ' paramvoigt.txt `


echo G height = $gheight
echo G location = $gloc
echo G width = $gwidth
echo G area = $garea

echo D1 height == area = $d1height
echo D1 location = $d1loc
echo D1 width = $d1width
echo D1 area = $d1area

echo D2 height == area = $d2height
echo D2 location = $d2loc
echo D2 width = $d2width
echo D2 area = $d2area

echo Lorentz solid offset = $bgint
echo Lorentz inclined offset = $bggrad

#Calculate SNR

snrmax=`awk 'BEGIN { max = -100000000000 } (($1 > 1700) && ($1 < 1800)) { if ( max < $2 ) max = $2 } END { print max }' bgremoved.xy`
snrmin=`awk 'BEGIN { min = 100000000000 } (($1 > 1700) && ($1 < 1800)) { if ( min > $2 ) min = $2 } END { print min }' bgremoved.xy`
bgremovedmax=`awk 'BEGIN { max = -100000000000 } ( $1 > 1200 ) { if ( max < $2 ) max = $2 } END { print max }' bgremoved.xy`

snr=`echo "scale=0; $bgremovedmax / ($snrmax - $snrmin)" | bc`
echo snrmax = $snrmax
echo snrmin = $snrmin
echo bgremovedmax = $bgremovedmax
echo snr = $snr

# R1 ratio
r1ratio=`echo "scale=5; $d1height / $gheight" | bc`
echo R1 ratio = $r1ratio

r2ratio=`echo "scale=5; $d1area / ( $d1area + $garea + $d2area )" | bc`
r2ratioa=`echo "scale=2; $d1area / ( $d1area + $garea + $d2area )" | bc`
r2ratiovoigt=$r2ratio
r2ratiocheck=`echo "scale=0; 100 * $r2ratioa / 1" | bc`
echo R2 ratio = $r2ratio
r2tempa=`echo "scale=3; ((-445 * $r2ratio) + 641 ) / 1 " | bc`
echo R2 temp = $r2tempa

ra2temp="na"
ra1temp="na"
d3loc=na
d3height=na
d3width=na
d3area=na
d4loc=na
d4height=na
d4width=na
d4area=na

fitstyle=Voigt2
temp="Standard R2 temp $r2tempa"
plottemp=$r2tempa

totalwidth=`echo "scale=2; $gwidth + $d1width + $d2width" | bc`

savefile=savefilevoigt.plt

else

fitstyle=Lorentzians
temp="Low temperature temp $ra2temp"
plottemp=$ra2temp

totalwidth=`echo "scale=2; $gwidth + $d1width + $d2width" | bc`

savefile=savefilelor.plt

fi

fi

## Double widths from HWHM to FWHM

totalwidth=`echo "scale=2; $totalwidth * 2" | bc`
gwidth=`echo "scale=2; $gwidth * 2" | bc`
d1width=`echo "scale=2; $d1width * 2" | bc`
d2width=`echo "scale=2; $d2width * 2" | bc`
d3width=`echo "scale=2; $d3width * 2" | bc`
d4width=`echo "scale=2; $d4width * 2" | bc`


echo $temp
echo $fitstyle

# Find the number of iterations
iterations="too many"
iterations=`awk ' /the fit converged/ { print $2 } ' fit.log`

# Go back to Gnuplot to make the final graphs

gnuplot $persist<<EOF

load "$savefile"

set title "$nicename R2 = $r2ratioa Temp = $plottemp"

# Uncomment the following to line up the axes
# set lmargin 6

set term x11 font ",6"

set size 1,1
set origin 0,0

set multiplot

set title "$nicename"

set size 0.4,0.5
set origin 0,0.5
plot [x=1000:1900] bg(x) title "Background", '$1' title "'$nicename'" with lines

set title "$nicename R2 = $r2ratioa Temp = $plottemp"

set size 0.6,1
set origin 0.4,0
plot [x=1000:1900] '$1' using 1:(\$2 - bg(\$1)) title "Background removed" with lines, d2peak(x), d1peak(x), gpeak(x), d3peak(x), d4peak(x), p(x) title "$fitstyle"
#plot [x=1000:1900] '$1' using 1:(\$2 - bg(\$1)) title "Background removed" with lines, d2peak(x), d1peak(x), gpeak(x), p(x)

set title "$nicename - $iterations iterations"

set size 0.4,0.5
set origin 0,0
plot [x=1000:1900] f(x) title "$fitstyle", '$1' title "'$nicename'" with lines, '$1' using 1:(\$2 - f(\$1)) title "Residual" with lines

unset multiplot
reset

set term post landscape color solid 8
set output 'combined.ps'

# Uncomment the following to line up the axes
# set lmargin 6

#set size ratio 1.5 1.5,1
set origin 0,0

set multiplot

set title "$nicename"

set size 0.33,0.5
set origin 0,0.5
plot [x=1000:1900] bg(x) title "Background", '$1' title "$nicename" with lines

set title "$nicename R2 = $r2ratioa RA2= $ra2ratio Temp = $plottemp"

set size 0.67,1
set origin 0.33,0
plot [x=1000:1900] '$1' using 1:(\$2 - bg(\$1)) title "Background removed" with lines, d2peak(x), d1peak(x), gpeak(x), d3peak(x), d4peak(x), p(x) title "$fitstyle"
#plot [x=1000:1900] '$1' using 1:(\$2 - bg(\$1)) title "Background removed" with lines, d2peak(x), d1peak(x), gpeak(x), p(x)

set title "$nicename after fitting - $iterations iterations"

set size 0.33,0.5
set origin 0,0
plot [x=1000:1900] f(x) title "$fitstyle", '$1' title "$nicename" with lines, '$1' using 1:(\$2 - f(\$1)) title "Residual" with lines

unset multiplot
reset

#set title "$nicename R2 = $r2ratioa 

set terminal png
set output 'data.png'
plot [x=1000:1900] bg(x) title "Background", '$1' title "'$nicename'" with lines

set output 'fit.png'
plot [x=1000:1900] f(x) title "$fitstyle", '$1' title "'$nicename'" with lines, '$1' using 1:(\$2 - f(\$1)) title "Residual" with lines

set output 'peaks.png'
plot [x=1000:1900] '$1' using 1:(\$2 - bg(\$1)) title "Background removed" with lines, d2peak(x), d1peak(x), gpeak(x), d3peak(x), d4peak(x), p(x) title "$fitstyle"
#plot [x=1000:1900] '$1' using 1:(\$2 - bg(\$1)) title "Background removed" with lines, d2peak(x), d1peak(x), gpeak(x), p(x)

EOF


# Output results to text files in case they're needed later
awk ' NR>4 && length > 3 { print $1, $2 / '$bgremovedmax' } ' lorentzians.xy > ${nicename}lorentzianschart.xy
awk ' NR>3 && length > 3 { print $1, $2 / '$bgremovedmax' } ' bgremoved.xy > ${nicename}bgremovedchart.xy

ps2pdf combined.ps ${nicename}combined.pdf
rm combined.ps

mv data.png ${nicename}data.png
mv fit.png ${nicename}fit.png
mv peaks.png ${nicename}peaks.png

#if (( "$snr" > "$threshold" ))
#then
#mv ${nicename}bgremovedchart.xy combinedfigure
#mv ${nicename}lorentzianschart.xy combinedfigure
#fi


#Output to combined results

echo $nicename $gheight $gloc $gwidth $garea $d1height $d1loc $d1width $d1area $d2height $d2loc $d2width $d2area $d3height $d3loc $d3width $d3area $d4height $d4loc $d4width $d4area $r1ratio $r2ratio $r2tempa $ra1ratio $ra1temp $ra2ratio $ra2temp $r2ratiovoigt $plottemp $totalwidth $totalwidthvoigt >> acombinedresults.txt



	echo Congratulations, new sample analysed 
	echo
return
else
	echo $outputresult
	echo Sample already processed
	echo
fi
}


function tidy {

#Tidy up

mv *.png jpg
mv *combined.pdf pdf
mv *chart.xy xy_chart_files
rm param*.txt
rm fit.log
rm *.xy
rm *.plt

#Make combined output
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=allanalysedraman.pdf pdf/*combined.pdf
}



while :
do
echo $# to go
if [[ "$#" > "0" ]]
then 

echo $# files left to process
processsample $1
shift

else 

tidy

exit
fi
done
