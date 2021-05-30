This is the code used to generate the results for the ICML 2021 conference paper "Fast Algorithms for Stackelberg Prediction Game with Least Squares Loss".

This code has the following requirements:
- MATLAB 2019b or later
- MOSEK solver (https://www.mosek.com/)
- The MATLAB Statistics and Machine Learning Toolbox (https://uk.mathworks.com/products/statistics.html)
- The MATLAB Parallel Computing Toolbox (https://uk.mathworks.com/products/parallel-computing.html)
- The MATLAB Optimization Toolbox (https://uk.mathworks.com/products/optimization.html)

All experimental results in our paper are saved as csv format in the file of "result". The file "datasets" is where our synthetic and real-world datasets are placed.  To save space, we only put wine and insurance datasets in "datasets" file, you can download other real-world datasets  from the following websites.

+ BlogFeedBack Dataset (https://archive.ics.uci.edu/ml/datasets/BlogFeedback)
+ Residential Building Dataset (https://archive.ics.uci.edu/ml/datasets/Residential+Building+Data+Set)

If you want to generate synthetic data, please run the ipynb file `make_regression_dataset.ipynb`. Then the correspondingly synthetic datasets would be placed at file "datasets/synthetic/".

### How to get results

+ To run the experiments of real datasets, please `run main_real_dataset_mosek.m`. 
+ To run the experiments of synthetic datasets, please `run main_synthetic_mosek.m`.

### How to plot the figures

+ To plot the average mean squared error (MSE) of different methods,  you can use the function `plot_mse(csvname)`, where "csvname" is the path of corresponding  MSE csv file. 
  For example, `csvname = './result/wine_modest_None_mse.csv'`
+ To obtain the figures of time comparison, you can use the function `plot_time(csvname)`, where "csvname" is the path of corresponding  time csv file. 
  For example, `csvname = './result/wine_modest_None_time.csv'`

### Citation

If you found the provided code useful, please cite our work.
```
@misc{wang2021fast,
      title={Fast Algorithms for Stackelberg Prediction Game with Least Squares Loss}, 
      author={Jiali Wang and He Chen and Rujun Jiang and Xudong Li and Zihao Li},
      year={2021},
      eprint={2105.05531},
      archivePrefix={arXiv},
      primaryClass={math.OC}
}
```
