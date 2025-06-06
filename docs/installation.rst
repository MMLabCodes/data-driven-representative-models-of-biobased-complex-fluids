Installation
============

This python package is designed to work within the unix shell. If you are working with **windows**, install miniconda on **ubuntu** https://ubuntu.com/ 

To get started you will need **anaconda** or **miniconda**:   

Anaconda: https://docs.anaconda.com/anaconda/install/   
Miniconda: https://docs.anaconda.com/miniconda/install/ (recommended for use with ubuntu)   

**Ubuntu** is a linux system emulator which is required to run some of the programmes utilised by this python package - please install miniconda using the ubuntu command line.

Creating an environment
-----------------------

To create an environment, open the command prompt (ubuntu if you are using windows) and execute the following:

.. code-block:: bash

	conda create --name DataDrivenModels
	conda activate DataDrivenModels

*Note: You do not have to call your environment *DataDrivenModels* as long as the name makes sense to you*

Now the packages and different dependancies can be installed. These are all located in the **environment.yml** file in the **docs** folder of the github repository.
This file can be used to install of the required pacakges.

.. code-block:: bash
	
	conda env update --file docs/environment.yml

Testing
-------

To ensure some key packages are installed it is important to run some tests. Execute the commands below and compare with the images.

**Antechamber** is part of AmberTools and carries out the parameterization of the molecules//polymers.

.. code-block:: bash
	
	Antechamber

If antechamber is available, you will see something similar to the following:

.. image:: images/antechamber.PNG

**Tleap** is an AmberTools programme that allows for the building of systems and generation of topology and coordinate files for 

.. code-block:: bash
	
	Tleap

If tleap is available, you will see something similar to the following:

.. image:: images/tleap.PNG

To exit tleap hold ctrl+c

A test for python and some associated packages is also important.

.. code-block:: bash
	
	python3

If python is available, you will see something similar to the following:

.. image:: images/python.PNG

*Note: this will open the python interpreter and you can enter any python commands after '>>>'*

Cloning the repository - Normal method
--------------------------------------

Now the repository needs to be cloned to give access to pyton modules within this pacakge. Go to the repository page in github https://github.com/MMLabCodes/data-driven-representative-models-of-biobased-complex-fluids and find the blue button labelled '<> code'.
Click here and select 'HTTPS' and copy the link. Now return to command line (ubuntu in windows), ensure you are in your home directory and execute the folowing:

.. code-block:: bash
	
	git clone copied link

This will clone the repository into your pc and you will be able to access all the required files. 
Don't forget you can navigate through the file explorer to view these files (if you are using windows, look for the linux folder with the penguin).

If this method doesn't work, see the alternative method below.


Cloning the repository - Alternative method
-------------------------------------------

First you will need to obtain a personal access token from github, once you have logged into github, click on your profile in the top right and navigate to (settings --> developer settings --> personal access tokens --> Tokens (classic)). 
Here, click on "generate new token --> generate new token (classic)" and enter a note "clone repo" and in the tick boxes, select "repo". 
Now scroll to the bottom and "generate token". This will give you a token you will need for the next step.

Now you can navigate to your home directory and execute the following commands:

.. code-block:: bash
	
	git clone https://USERNAME:YOUR_TOKEN@github.com/MMLabCodes/data-driven-representative-models-of-biobased-complex-fluids.git
	cd data-driven-representative-models-of-biobased-complex-fluids

The final 'cd' command navigate to the directory containing the notebooks and scripts required for the tutorials.

Jupyter Notebook tutorials
==========================

There are a couple of jupyter notebooks that contain tutorials and are aptly labelled. A section can be found explaining the contents of each notebook.
To launch jupyter notebooks execute the following in the command line (ubuntu in windows):

.. code-block:: bash
	
	jupyter notebook

This will launch a local jupyter notebook server and a series of URLs will be returned. Copy the **first link** containing ('localhost:8888') and copy//paste into a browser.
From there, refer to the section in this documentation about different tutorials.

Even though there are pre-existing jupyter notebook notebook tutorials, there is enough guidance and explantion in the following section of this documentation
to carry out all the tasks within this python package related to the generation of data-driven-representative-models-of-biobased-complex-fluids.

1_Generating_models_of_complex_fluids
-------------------------------------
The contents of this notebook is explained thoroughly in the **Guide to Model Generation** section of this documentation.
It covers:

- How to structure raw data for use with this package
- How to generate representative models of complex fluids
- How to evaluate the efficacy of all generated models

This notebook has good explanations throughout and is the place to to start with using this code.

1b_Generating_models_of_complex_fluids_QUICKSTART
-------------------------------------------------
This notebook contains the exact some code as **1_Generating_models_of_complex_fluids** but with entirely stripped back explanations.

This is the place to start if you are very experienced with python or would like to adapt this code into python scripts for use with command line.


2_Complex_fluid_models_to_amber
-------------------------------

EXPLANATION AND GUIDE COMING SOON...

.. note:: 

   Whilst no guide currently exists in the documentation this jupyter notebook contains enough information to guide you through the process of generating
   amber files for complex fluid model systems.


2b_Complex_fluid_models_to_amber_QUICKSTART
-------------------------------------------

EXPLANATION AND GUIDE COMING SOON...










