java -jar ../../EVA.jar ./BEST_01/configuration.txt r=10 technique=GB w=tumor_weights_RFCI.txt
mv ../../data ./BEST_01/data

java -jar ../../EVA.jar ./BEST_04/configuration.txt r=10 technique=GB w=tumor_weights_RFCI.txt
mv ../../data ./BEST_04/data

java -jar ../../EVA.jar ./BEST_07/configuration.txt r=10 technique=GB w=tumor_weights_RFCI.txt
mv ../../data ./BEST_07/data


java -jar ../../EVA.jar ./WORST_01/configuration.txt r=10 technique=GB w=tumor_weights_RFCI.txt
mv ../../data ./WORST_01/data

java -jar ../../EVA.jar ./WORST_04/configuration.txt r=10 technique=GB w=tumor_weights_RFCI.txt
mv ../../data ./WORST_04/data

java -jar ../../EVA.jar ./WORST_07/configuration.txt r=10 technique=GB w=tumor_weights_RFCI.txt
mv ../../data ./WORST_07/data