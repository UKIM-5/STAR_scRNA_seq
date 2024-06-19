# STAR_Scheiber_et_al_2024

The following script allows you to reproduce the results and plots from our paper "Profiling low-mRNA content cells 
in complex human tissues using BD Rhapsody single-cell analysis".
For the sake of simplicity, we hereby only include basic processing steps such as data read-in, 
quality control as well as manual cell type annotation on characteristic markers. 
The script can be run via the terminal for all three sample types (liver, lung and prostate).
Original files and quality thresholds are taken from our previously published studies (Salcher et al 2022, Hautz et al 2023, Heidegger et al 2022), 
where we showed the successful identification of low mRNA cells via scRNA-seq on bigger datasets.
We highly recommend checking out the official tutorial of scanpy (https://scanpy.readthedocs.io) as well as the single-cell best practices continously 
updated by the Theis group (https://www.sc-best-practices.org).

As this approach includes slightly more steps and code, we also included a jupyter notebook as hmtl which only uses the code available
in the protocol itself (STAR_Protocols_Code_PAPER). The raw files are available on Zenodo including the html file with the code of the paper.

NOTE: The following setup is only usable with a virtual conda environment, thus needing Anaconda/Minoconda installment beforehand.
Please refer to the corresponding website for more information: https://www.anaconda.com

For an easy setup, please read and follow the next steps:

1. Please download or clone this repository.

2. Navigate to the your current directory via the terminal and install a new virtual conda environment using the following command:
$ conda env create -f environment.yml

3. Activate the environment using:
$ conda activate STAR_env

4. Navigate to the raw_data directory and get the h5ad files from zenodo: (10.5281/zenodo.12165000)

5. Upon raw data is deposited, go to the parent folder where the main_function is deposited and call 
the script - choose which sample you want to use the script for (lung, prostate or liver):
$ python star_main_function.py --sample_type=lung

NOTE: Scanpy requires > python 3.6, however to our knowledge most higher python versions are working fine (we used 3.12 in our case). 

Output:
Adata file including all figures used for cell annotation as well as final plots presented in Scheiber et al 2024.
