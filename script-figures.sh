#!/bin/sh

# @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
# Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
# @author Bertrand Kerautret (\c kerautre@liris.cnrs.fr )
# Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France

# Define the path to the program
EXEC=build/bin/tv-triangulation-color

# Choose your zoom factor
ZOOM=16

# Choose your image input and output path
INPUTPATH=Input
OUTPUTPATH=Output

###############################################################################
echo "------ Creating Figure 1 ----------------------------------------------"
IMAGES="x-shapes.pgm dolphin.ppm ipol-coarsened.pgm barbara.ppm ara.ppm"
let n=1
for i in ${IMAGES}; do
    INPUT=${INPUTPATH}/$i
    FILENAME=`basename $i`
    OUTPUT="${OUTPUTPATH}/${FILENAME%.*}-x${ZOOM}"
    echo "- Creating Fig. 1, col. 1, row $n: ${OUTPUT}-bitmap.png"
    convert -filter Point -resample "${ZOOM}00%" ${INPUT} ${OUTPUT}-bitmap.png
    echo "- Creating Fig. 1, col. 2, row $n: ${OUTPUT}-lg-nogtv-smooth.png"
    ${EXEC} -i ${INPUT} -b ${ZOOM} -S 0 -L 0 -o "tmp" 1> tmp.log 2> tmp.err
    mv "tmp-lg.png" "${OUTPUT}-lg-nogtv-smooth.png"
    echo "- Creating Fig. 1, col. 3, row $n: ${OUTPUT}-lg-smooth.png"
    ${EXEC} -i ${INPUT} -b ${ZOOM} -S 0 -o "tmp"  1> tmp.log 2> tmp.err
    mv "tmp-lg.png" "${OUTPUT}-lg-smooth.png"
    echo "- Creating Fig. 1, col. 4, row $n: ${OUTPUT}-lg-crisp.png"
    ${EXEC} -i ${INPUT} -b ${ZOOM} -S 100 -o "tmp"  1> tmp.log 2> tmp.err
    mv "tmp-lg.png" "${OUTPUT}-lg-crisp.png"
    let n=n+1
done

###############################################################################
echo "------ Creating Figure 2 ----------------------------------------------"
IMAGES="dolphin.ppm"
let n=1
for i in ${IMAGES}; do
    INPUT=${INPUTPATH}/$i
    FILENAME=`basename $i`
    OUTPUT="${OUTPUTPATH}/${FILENAME%.*}-dual"
    echo "- Creating Fig. 2, col. 1, row 1: ${OUTPUT}-nogtv.eps"
    ${EXEC} -i ${INPUT} -D 0 -o "tmp" -L 0 -R 0 -E ${OUTPUT}-nogtv.eps  1> tmp.log 2> tmp.err
    echo "- Creating Fig. 2, col. 2, row 1: ${OUTPUT}-nogtv-reg.eps"
    ${EXEC} -i ${INPUT} -D 0 -o "tmp" -L 0 -E ${OUTPUT}-nogtv-reg.eps  1> tmp.log 2> tmp.err
    echo "- Creating Fig. 2, col. 1, row 2: ${OUTPUT}-noreg.eps"
    ${EXEC} -i ${INPUT} -D 0 -o "tmp" -R 0 -E ${OUTPUT}-noreg.eps  1> tmp.log 2> tmp.err
    echo "- Creating Fig. 2, col. 2, row 2: ${OUTPUT}.eps"
    ${EXEC} -i ${INPUT} -D 0 -o "tmp" -E ${OUTPUT}.eps  1> tmp.log 2> tmp.err
    let n=n+1
done

###############################################################################
echo "------ Creating Figure 3 ----------------------------------------------"
IMAGES="dolphin.ppm"
let n=1
for i in ${IMAGES}; do
    INPUT=${INPUTPATH}/$i
    FILENAME=`basename $i`
    OUTPUT="${OUTPUTPATH}/${FILENAME%.*}-x${ZOOM}"
    echo "- Creating Fig. 3, col. 1: ${OUTPUT}-2nd-with-lines-crisp.png"
    ${EXEC} -i ${INPUT} -b ${ZOOM} -D 32 -S 100 -o "tmp" 1> tmp.log 2> tmp.err
    mv "tmp-2nd-with-lines.png" "${OUTPUT}-2nd-with-lines-crisp.png"
    echo "- Creating Fig. 3, col. 2: ${OUTPUT}-2nd-crisp.png"
    ${EXEC} -i ${INPUT} -b ${ZOOM} -D 16 -S 100 -o "tmp" 1> tmp.log 2> tmp.err
    mv "tmp-2nd.png" "${OUTPUT}-2nd-crisp.png"
    let n=n+1
done

###############################################################################
echo "------ Creating Figure 4 ----------------------------------------------"
IMAGES="x-shapes.pgm dolphin.ppm ipol-coarsened.pgm barbara.ppm ara.ppm"
let n=1
for i in ${IMAGES}; do
    INPUT=${INPUTPATH}/$i
    FILENAME=`basename $i`
    OUTPUT="${OUTPUTPATH}/${FILENAME%.*}-x${ZOOM}"
    echo "- Creating Fig. 4, col. 4, row $n: ${OUTPUT}-2nd.png"
    time ${EXEC} -i ${INPUT} -b ${ZOOM} -D 16 -o "tmp" 1> tmp.log 2> tmp.err
    mv "tmp-2nd.png" "${OUTPUT}-2nd.png"
    let n=n+1
done

# Cleaning stuff
rm tmp.log tmp.err output-laplacian-before.pgm output-tv.ppm

