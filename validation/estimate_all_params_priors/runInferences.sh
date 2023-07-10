module load openjdk/17.0.0_35
for seed in {101..200}; do
cd "run_${seed}/inf";
bsub -J "prior${seed}" -W 24:00 -R "rusage[mem=500]" -n 4 java -Dglass.platform=Monocle -Dmonocle.platform=Headless --module-path ../../../../javafx-sdk-17.0.6-linux-monocle/lib --add-modules javafx.fxml,javafx.base -jar ../../../../stratigraphic-ranges.jar -seed ${seed} -threads -1 -version_file ../../../../version-ranges.xml -version_file ../../../../version-feast.xml -version_file ../../../../version-beast.xml -version_file ../../../../version-sa.xml -version_file ../../../../version-beastLabs.xml -version_file ../../../../version-mm.xml -version_file ../../../../version-beastfx.xml  -overwrite sRanges_inference.xml
cd ../../;
done
