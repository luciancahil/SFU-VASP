CRYSTAL=$1
SETTINGS=$2
timestamp=$(date +%s)
mainFolder=$CRYSTAL$timestamp


if [ -z "$CRYSTAL" ]; then
  echo "Variable CRYSTAL is empty or unset"
  exit 0
fi
if [ -z "$SETTINGS" ]; then
  echo "Variable SETTINGS is empty or unset"
  exit 0

fi

mkdir jobs/$mainFolder -p
echo $mainFolder

cd jobs
cd $mainFolder
cp ../../vasp.sh vasp.sh

echo cd  jobs/$mainFolder 
echo bash vasp.sh $CRYSTAL $SETTINGS
sbatch vasp.sh $CRYSTAL $SETTINGS
