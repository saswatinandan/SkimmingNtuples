#!/bin/sh
if [ ! $1 ] ; then
echo "Please specify the code you want to compile by typing :"
echo "./Make <Your-Code.C>"
exit 1
fi
echo "================================================================"

echo "====> Here are the libraries:" "\n" `root-config --cflags --glibs` -lTMVA "\n"
filename=`echo $1 | awk -F"." '{print $1}'`
exefilename=${filename}.exe
rm -f $exefilename
#g++ $1  -o $exefilename `root-config --cflags --glibs` -lTMVA
g++ $1 -I$CMSSW_BASE/src `root-config --cflags` -o $exefilename -L$CMSSW_BASE/lib/$SCRAM_ARCH/ -lHTT-utilitiesLepEffInterface `root-config --libs` -lMinuit -lTMVA -lm
echo ""
if [ -e $exefilename ]; then
echo "====> Created exe file : "
ls -lrt $exefilename
echo "====> Done."
else
echo "====> Did not create the exe file!"
fi
echo "================================================================"

