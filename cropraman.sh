while :
do
echo $# to go

if [[ "$#" > "0" ]]
then 
echo "cropping file $1"
`awk ' ( $1 > 800 ) && ( $1 < 2200 ) { print $1, $2 } ' $1 > crop$1`

shift


else
echo finish
exit
fi
done
