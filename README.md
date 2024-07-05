# Polish

## Przygotowanie środowiska

Poniższe kroki umożliwiają przygotowanie środowiska do przeprowadzenia analizy WGCNA:
git clone https://github.com/SleepDealler/WGCNA.git
cd ./WGCNA/
chmod +x libs.R
./libs.R

## Dane

Dane wymagane do przeprowadzenia analizy:
- Transponowana macierz zliczeń (nazwy kolumn jako geny, nazwy wierszy jako nazwa próbki). Wywoływane  z argumentem -i.
- Metadane (wartości w kolumnie 'Library Name' identyczne do nazwy wierszy macierzy zliczeń. Wywoływane  z argumentem -m.
- Plik z rozszerzeniem csv, w którym kolejno wymienionymi wartościami są cechy filogenetyczne. Wywoływane z argumentem -p.
Przykładowe dane dostępne w ./example.

## Argumenty przyjmowane przez skrypt

Wymagane argumenty:
-i lub --input (string): Ścieżka do pliku wejściowego zawierającego dane ekspresji genów (CSV).
-m lub --metadata (string): Ścieżka do pliku metadanych zawierającego informacje o próbkach (CSV).
-p lub --phylo (string): Ścieżka do pliku z cechami filogenetycznymi (CSV).
Opcjonalne argumenty:
-d lub --deepsplit (integer): Głębokość podziału dla dynamicznego cięcia drzewa. Domyślna wartość to 2.
-mc lub --minClustersize (integer): Minimalna wielkość klastra. Domyślna wartość to 20.
-ch lub --cutHeight (double): Wysokość cięcia drzewa. Domyślna wartość to 0.99.
-mdt lub --MeDissThres (double): Próg podobieństwa dla scalania modułów eigengenów. Domyślna wartość to 0.25.

Przykładowe użycie:
./WGCNA.R -i ./example/input.csv -m ./example/metadata.csv -p ./example/phylo.csv

## Wyniki

Skrypt zwraca wykresy z poszczególnych etapów analizy (./plots/) oraz pliki niezbędne do konstrukcji sieci w Cytoscape (./results/)

## WGCNA_v2.Rmd

W przypadku braku możliwości wykonania analizy przez WGCNA.R, możliwe jest przeprowadzenie analizy przez WGCNA_v2.Rmd na domyślnych parametrach. Należy w skrypcie zmienić nazwy otwieranych plików lub dostosować nazwy plików do skryptu (input.csv, metadata.csv, phylo.csv). Podczas wykonywania analizy przez WGCNA_v2.Rmd wcześniej powinien być uruchomiony Cytoscape.

# English

## Environment Setup

The following steps allow you to set up the environment to perform WGCNA analysis:
git clone https://github.com/SleepDealler/WGCNA.git
cd ./WGCNA/
chmod +x libs.R
./libs.R

## Data

Data required to perform the analysis:
- Transposed count matrix (column names as genes, row names as sample names). Called with the argument -i.
- Metadata (values in the 'Library Name' column identical to the row names of the count matrix). Called with the argument -m.
- A CSV file listing phylogenetic traits in sequence. Called with the argument -p.
Example data available in ./example.

## Script Arguments

Required Arguments:
-i or --input (string): Path to the input file containing gene expression data (CSV).
-m or --metadata (string): Path to the metadata file containing sample information (CSV).
-p or --phylo (string): Path to the CSV file containing phylogenetic traits.
Optional Arguments:
-d or --deepsplit (integer): Depth of split for dynamic tree cutting. Default value is 2.
-mc or --minClustersize (integer): Minimum cluster size. Default value is 20.
-ch or --cutHeight (double): Tree cut height. Default value is 0.99.
-mdt or --MeDissThres (double): Module eigengene dissimilarity threshold for merging modules. Default value is 0.25.

Example usage:
./WGCNA.R -i ./example/input.csv -m ./example/metadata.csv -p ./example/phylo.csv

## Results

The script outputs plots from various stages of the analysis (./plots/) and files required to construct networks in Cytoscape (./results/).

## WGCNA_v2.Rmd

If the analysis cannot be performed using WGCNA.R, it can be conducted using WGCNA_v2.Rmd with default parameters. You need to change the file names in the script or adapt the file names to the script (input.csv, metadata.csv, phylo.csv). Ensure Cytoscape is running before performing the analysis with WGCNA_v2.Rmd.