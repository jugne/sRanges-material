module load openjdk/17.0.0_35
for seed in {1..100}; do
cd "run_${seed}/inf";
bsub -J "sr_${seed}" -W 24:00 -R "rusage[mem=10000]" java -Dglass.platform=Monocle -Dmonocle.platform=Headless -Dprism.order=sw  --module-path ../../javafx-sdk-17.0.6-linux-monocle/lib --add-modules javafx.fxml,javafx.base -jar ../../stratigraphic-ranges.jar -seed ${seed} -version_file ../../version-ranges.xml -version_file ../../version-feast.xml -version_file ../../version-beast.xml -version_file ../../version-sa.xml -version_file ../../version-beastLabs.xml -version_file ../../version-mm.xml -version_file ../../version-beastfx.xml  -overwrite sRanges_inference.xml
cd ../../;
done
