
Male pedigree toolbox
=====================

This is a collection of functionalities for the analysis of male pedigrees based on Y-chromosomal markers. Here follows
a short overview of functionalities:


* Generational distance calculations between all individuals in a pedigree.
* Number of mutations between all alleles for all markers in a pedigree.
* Infer alleles and mutation events in pedigrees and draw these pedigrees.
* Cluster alleles/individuals based on mutation distance between them
* Simulate mutations based on marker mutation rates and use these simulations to train various machine learning models for the prediction of generational distance between individuals based on markers.

Contents
--------


* Installing

  * `Executables <#download-executable>`_
  * `Python package installation <#clone-and-pip-install>`_
  * `Execute from main script <#execute-from-main>`_

* Running

  * Pedigree investigatioN

    * `Meiotic distances <#meiotic-distances-in-pedigrees-distance>`_
    * `Counting mutations between alleles <#counting-mutations-between-alleles-of-markers-mut_diff>`_
    * `Inferring pedigree mutation events <#infering-pedigree-mutation-events-ped_mut_graph>`_
    * `Clustering alleles <#clustering-alleles-based-on-mutation-distance-draw_pedigrees>`_
    * `Run it all <#run-all-the-above-commands-in-tandem-all>`_

  * Generational distance prediction

    * `Simulate data <#simulate-alleles-data-simulate-command-line-only>`_
    * `Make models <#create-classification-models-from-simulated-data-make_models-command-line-only>`_
    * `Predict generational distance <#predict-generational-distance-predict>`_

* `Full example <#full-example>`_

Installing:
===========

Download executable
-------------------

The easiest way of using the Male pedigree toolbox is by using the precompiled executables that have been created for
linux and windows. Unfortunately there is no executable for mac available. The downside of these executables is that it
takes a long time for them to start up (around 20 seconds). There is a gui and command line executable available.

In order for inferring pedigree mutations to work properly https://graphviz.org/ is required. You will need to add Graphiz
to your system path or add the Graphiz directory inside of the MalePedigreeToolbox folder. Graphiz is included in the executables.

Clone and pip install
^^^^^^^^^^^^^^^^^^^^^

The repository can also be installed with pip for convenient command line acces. This also allows you to start the
graphical user interface from the command line. In order for the tool to be able to start python 3.6 or higher is
required

Installing with pip is as simple as :

.. code-block::

   $ pip install male-pedigree-toolbox

This will install this toolbox as a python package and make it available on the command line. Now check that the
command line interface of the toolbox is properly installed:

.. code-block::
   # print the current version, this is an example the number will differ
   $ mpt --version
   MalePedigreeToolbox 0.1.2

You can check the same for the GUI. This command should start up a GUI.

.. code-block::

   $ mpt_gui

Execute from main
^^^^^^^^^^^^^^^^^

In case the executable does not work, and you don't want to pip-install the package. You can always clone the GitHub
repository and execute the main.py script:

.. code-block::

   $ clone https://github.com/genid/MalePedigreeToolbox.git
   $ python main.py --version
   MalePedigreeToolbox v0.1.0-beta

Or navigate into the `gui <./MalePedigreeToolbox/gui>`_ folder and execute the main_gui.py script:

.. code-block::

   $ python main_gui.py

Keep in mind that the following python packages are required as well as python 3.6:


* pandas
* numpy
* statsmodels
* scipy
* matplotlib
* joblib
* sklearn
* scikit-learn
* tqdm
* openpyxl
* graphviz

This package is required for the gui:


* PySimpleGUI

All of these packages can be installed with pip:

.. code-block::

   # one at a time
   $ pip install <package name>
   # or all at once
   $ pip install -r requirements.txt

Running
=======

There are a number of different functionalities that can be used from this toolkit. Here follows an explanation for each
of these functionalities with some example in and outputs. The examples are for the command line but the same applies
for the inputs of the GUI unless stated otherwise. Alternatively you can always make use of -h or --help to get an
overview of all options available for a certain subcommand.

Pedigree investigation commands
-------------------------------

These are commands that can be used to investigate pedigrees in a number of ways.

Meiotic distances in pedigrees (distance)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Calculate distances between all individuals in the provided pedigrees. The pedigrees need to be in Trivial
Graph Format (tgf). The command can calculate the distances between all individuals in a pedigree.

Example command:

.. code-block::

   $ mpt distances -t example_tgfs -o output_folder

This will create a comma separated values (csv) file containing the generational distance between all individuals of
each pedigree in the specified output folder.

Counting mutations between alleles of markers (mut_diff)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Get the number of mutations between all alleles for all markers in pedigrees. The input for this command is an alleles
file. This is a .csv file that contains the alleles for each marker of one or more pedigrees. An Example of an alleles
filecan be found at `examples/example_alleles.csv <./examples/example_alleles.csv>`_. The number of alleles does not have
to be 6. Optionally the distances between all individuals of the different pedigrees can be provided
(this can be generated with the `distance <#meiotic-distances-in-pedigrees-distance>`_ command).

Example command:

.. code-block::

   $ mpt mut_diff -af example_alleles.csv -df distances_output_folder/distances.csv -o output_folder

This always results in at least 2 files. Firstly, a full output file containing the number of mutations that occured
between all individuals of a pedigree for all markers for each allele. Secondly, a summary output file that takes the mutations for
all markers together and shows the number of mutations between all individuals of a pedigree. If a distance file was
specified then percentage of mutation is calculated for each number of meiosis present in the provided pedigrees. The -pf
flag can be specified as well to generate a file that can be used to simulate data for creating machine learning models
for the prediction of generational distance.

Infering pedigree mutation events (ped_mut_graph)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Infer alleles and mutation events for pedigrees containing individuals with unknown alleles. The input for this command
is an alleles file (for an example see the `mut_diff <#counting-mutations-between-alleles-of-markers-mut_diff>`_
description) and a folder containing pedigrees in .tgf format.

Example command:

.. code-block::

   $ mpt ped_mut_graph -af example_allelles.csv -t example_tgfs -o output_folder

This will generate a pedigree for each marker containing the number of mutations that occured between descendants in the
pedigree. It will also contain an overview graph for each pedigreewhere all unique sets of alleles get their own color.
Each pedigree also gets a file with mutation rates for each marker based on that pedigree. Finally, a file that summarizes
all these mutation rates for all pedigrees is also generated.


.. image:: ./examples/marker_example.png
   :target: ./examples/marker_example.png
   :alt: plot

*Example of a pedigree for a certain marker with inferred mutation locations. The number at the edge indicates the number
of mutations the color indicates where this mutation could have occured, since these mutations are annotated at the
first place that they could have occured.*


.. image:: ./examples/all_marker_example.png
   :target: ./examples/all_marker_example.png
   :alt: plot

*Example of the same pedigree for all markers. Here Each unique allele gets a unique color. A .csv file acompanies this
file giving information on what marker mutated on what edge. All edges where mutations occured have an id together with
the number of mutations that occured. Keep in mind that these mutations are placed at the first edge they
could have occured.*

Clustering alleles based on mutation distance (draw_pedigrees)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Identify likely related individuals based on the mutation distance of the alleles of measured markers. The input for
this functionality is full list of mutation distances between all markers for all alleles (this can be generated with
the `mut_diff <#counting-mutations-between-alleles-of-markers-mut_diff>`_ command). . For an example of
a mutation rates file see `examples/example_marker_rates.csv <examples/example_marker_rates.csv>`_. Additionally, for more
accurate results you can also provide the mutation rates for all markers in a separate file. You can either define the
number of clusters yourself or let the program calculate the optimal number using silhouette score to measure how
good the clustering is.

Example command:

.. code-block::

   $ mpt draw_pedigrees -fm mut_diff_out/full_out.csv -mr example_marker_rates.csv -o output_folder -t both

This will produce a dendrogram or multi-dimensional scaling (MDS) plot or both for each pedigree present in the full
mutation distances file. Besides that text files are provided that contain the clusters, in order to easily work with
get all the individuals of a certain cluster.

Run all the above commands in tandem (all)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There is a command to run all the above functionalities in order where files created from one command are used as inputs
for others. This requires at the minimum a folder with .tgf files and an alleles file to run.

Example command:

.. code-block::

   $ mpt all -af allele_file.csv -t tgf_folder -o output_folder

Pedigree prediction functions
-----------------------------

These are a set of commands that can be used to generate models for the prediction of generational difference between
based on the number of mutations one individual has compared to another.

Simulate alleles data (simulate) (command line only)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Simulate data for creating classification models based on mutation rates of markers. These mutation rates can be
obtained from `ped_mut_graph <#infering-pedigree-mutation-events-ped_mut_graph>`_ or calculated yourself. For an example of
a mutation rates file see `examples/example_marker_rates.csv <examples/example_marker_rates.csv>`_. This command
generates data for the `make_models <#create-classification-models-from-simulated-data-make_models-command-line-only>`_
command in order to have a sufficiently large dataset to create the models from. You can specify the number of
generations and the number of inidividuals per generation that you want to simulate. Each generation is simulated
independant from previous generations.

Example command:

.. code-block::

   $ mpt simulate -i marker_rate_file.csv -o output_folder -n 10000 -g 50

This will generate one file containing the simulated mutations for each marker of each individual
over all generations. We recommend generating for at least 10.000 individuals per generation. An example of  the
simulated data can be found at `examples/example_simulated.csv <./examples/example_simulated.csv>`__.

Create classification models from simulated data (make_models) (command line only)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create classification models that predict a generational distance between 2 individuals of 1 till the number of
simulated generations. There are a number of different models that can be chosen from. From our experience the best
performing models are the multi-layer perceptron, support vector machines (SVM, scale very badly with large datasets) and
linear discriminant analysis (LDA). Depending on the model this can run for quite a while. It is also advised to use a
large number of cores if available to speed up the calculations.

Example command (this command runs for a long time):

.. code-block::

   $ mpt make_models -i simulated_data.csv -o output_folder -mt MLP LDA -c -1

This will create a pickled RandomizedSearchCV object containing the model. These can be used by the final component of
these comands to predict the generational distance between individuals.

Predict generational distance (predict)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Allows to predict the generational distance between one or more individuals based on the number of mutations between a
sets of markers. There are a number of pre-computed models that can be used for a few standard sets of markers. The
following marker sets have pre-computed models:


* RMPLEX
* PPY23
* YFP
* PPY23 + RMPLEX
* YFP + RMPLEX

The input
file can be generated from an alleles file with the help of the
`mut_diff <#counting-mutations-between-alleles-of-markers-mut_diff>`_ command. The file should look the same as the
`examples/example_simulated.csv <./examples/example_predict_input.csv>`__.

Example command:

.. code-block::

   $ mpt predict -i marker_mutation_observations.csv -o output_folder -m model_file.joblib -tf simulated_data.csv

Full example
============

Here is an example for using the `all <#run-all-the-above-commands-in-tandem-all>`_ command using files provided in the
`examples <./examples>`_ folder of this repository. The example is for the command line specifically but the provided
output should be the same for the gui. Take note that the example command assumes that it is executed from MalePedigreeToolbox
base folder.

.. code-block::

   $ mpt all --tgf_folder ./examples/example_tgfs/ --allele_file ./examples/example_alleles.csv --outdir ./output_directory --type both --random_state 5 --marker_rates ./examples/example_marker_rates.csv --clusters 2

    INFO 15:11:57.464672 (0.004 sec) - Loading libraries...
    INFO 15:12:04.859765 (7.399 sec) - Running all modules in tandem...
    INFO 15:12:04.859927 (7.399 sec) -
    INFO 15:12:04.859969 (7.399 sec) - Step 1/4
    INFO 15:12:04.860012 (7.399 sec) - Started with calculating pairwise distances.
    INFO 15:12:04.861764 (7.401 sec) - Finished calculating pairwise distances
    INFO 15:12:04.861858 (7.401 sec) -
    INFO 15:12:04.861897 (7.401 sec) - Step 2/4
    INFO 15:12:04.861940 (7.401 sec) - Starting with calculating differentiation rates
    INFO 15:12:04.870831 (7.410 sec) - Finished reading both input files
    INFO 15:12:04.871125 (7.411 sec) - In total there are 49 markers that will be analysed.
    WARNING 15:12:04.872397 (7.412 sec) - Marker (DYS1001) is not present in 1036648 and 1992767. The comparisson will be skipped.
    INFO 15:12:04.874018 (7.413 sec) - Calculation progress: 45%...
    INFO 15:12:05.259639 (7.799 sec) - Starting with writing mutation differentiation information to files
    INFO 15:12:05.311114 (7.851 sec) - Started with summarising and writing meiosis differentiation rates to file
    INFO 15:12:05.323806 (7.863 sec) - Finished calculating differentiation rates.
    INFO 15:12:05.328681 (7.868 sec) -
    INFO 15:12:05.328735 (7.868 sec) - Step 3/4
    INFO 15:12:05.328782 (7.868 sec) - Starting with creating dendograms based on mutation differentiation
    INFO 15:12:05.841020 (8.380 sec) - Calculation progress: 100%...
    INFO 15:12:05.841089 (8.380 sec) - Finished drawing dendograms for all pedigrees that were present
    INFO 15:12:05.841177 (8.381 sec) -
    INFO 15:12:05.841211 (8.381 sec) - Step 4/4
    INFO 15:12:05.841253 (8.381 sec) - Start with caclulating mutations from pedigrees
    INFO 15:12:05.843494 (8.383 sec) - Processing pedigree 1
    INFO 15:12:09.122974 (11.662 sec) - Processing pedigree 73
    INFO 15:12:13.255955 (15.795 sec) - Calculation progress: 45%...
    INFO 15:12:13.276336 (15.816 sec) - Finished calculating mutations from pedigrees
    INFO 15:12:13.276672 (15.816 sec) - Finished running all modules
    INFO 15:12:13.276811 (15.816 sec) - The log file can be found at './run.log'

This will create all the files in a folder called output_directory located in the folder from which this command was
executed as well as a run.log file containing similar information to what was put on the command line.

Some additional notes when running analysis with the example files, when running the classifier_predict command with
the example file you should use the YFP model since the example data is for that model.
