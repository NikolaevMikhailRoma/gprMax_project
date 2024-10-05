GPRMax Project: Neural Networks for Ground Penetrating Radar Data Processing
This project implements methods for solving inverse problems in ground penetrating radar (GPR) data processing using neural networks.
Project Overview
The main objective is to process radargrams into electro-physical cross-sections. We focus on developing a simple, automated approach to GPR data interpretation.
Dataset Generation
Due to the difficulty of obtaining large volumes of labeled GPR data, we utilize the gprMax package on GPU to generate synthetic datasets. The simulation environment is set up on a local machine with PyCharm, CUDA, and a virtual environment. Two datasets have been created:

A simple dataset
A more complex dataset with added heterogeneities

Data Processing
The dataset was processed to address dimensionality issues (126x125 instead of 128x128). Additional features were extracted using the librosa module.
Neural Network Models
We have implemented and compared several standard image processing models, including:

U-NET
Autoencoder (AE)
PSP-Net
Simplified U-NET (based on a published architecture)
Seq_NET
Convolutional Autoencoder (CAE)
Variational Convolutional Autoencoder (VCAE)

Project Structure

ML.ipynb: Main notebook for neural network training and evaluation
Dataset/: Simple 2D model dataset
Dataset_2/: Complex 2D model dataset
Dataset_3/: Additional generated dataset
gprmax_library.py: Library of functions and classes for forward problem modeling
Tests/test_gprmax_library.py: Unit tests for the modeling library
generation_models.py: Code for model generation
gprMax_colab_template.ipynb: Template for running gprMax in Google Colab
terminal.txt: Quick reference for terminal commands

Implemented Features

Inverse problem solution
Horizontal resolution enhancement
Deconvolution (currently non-functional)

Known Issues and Future Work

Incorporate conductivity information into the CAE model
Adapt VAE architecture for conductor detection
Implement genetic algorithms for the simplified U-NET architecture
Address issues with Dataset_2 preprocessing
Develop LSTM networks for simple datasets with varying time dimensions

Data Availability
The datasets are stored in the cloud due to their large size.
Technical Details

The project uses various deep learning architectures and signal processing techniques.
GPU acceleration is employed for both data generation and model training.
Custom data preprocessing steps are implemented to handle dimensionality issues.

Conclusion
This project demonstrates the potential of neural networks in automating GPR data interpretation, achieving results that surpass manual processing capabilities.
For more detailed information, please refer to the presentation and results in Guide/ML_gpr.pptx.


