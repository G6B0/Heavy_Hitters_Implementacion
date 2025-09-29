# Heavy_Hitters_Implementacion
Diseño e implementación para la detección de heavy hitters de k-mers en un conjunto de genomas de Escherichia Coli utilizando dos estructuras de CountSketch (CS) y TowerSketch con CountMin y Conservative Update (TS-CMCU)
---------------------------------------------------
Para poder ejecutar los archivos usados en este repositorio se deben seguir las siguientes instrucciones:

Ejecutar act1
1) Clonar el repositorio, luego abrir una terminal en la ubicación /Heavy_Hitters_Implementacion
2) Compilar act1 : g++ act1.cpp -o act1
3) Darle permiso al .sh con el comando chmod +x ./scrip_act1.sh
4) Ejecutar scrip_act1.sh de la siguiente manera: ./scrip_act1.sh
5) En caso de no poder ejectuarlo debe instalar dos2unix (sudo apt-get install dos2unix)
6) Ejecutar dos2unix scrip_act1.sh y volver a paso 4
------------------------------------------------------
Ejecutar act3
1) Compilar act3 : g++ act3.cpp TowerSketch.cpp CountMin_CU.cpp CountSketch.cpp -o act3
2) Darle permiso al .sh con el comando chmod +x ./scrip_act3.sh
3) Ejecutar scrip_act3.sh
