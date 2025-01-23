### Python Tests

To install pytest:
````
pip install -U pytest
````

#### Run tests:  
Add the ``--not-cleanup`` argument to ensure that test directories and test files are not removed and can be used for debugging purposes.
   
All the provided arguments below can be used together to run all tests in 1 pytest run:

##### Test full _md_run_ pipline
````
pytest  Miniconda3/envs/md/lib/python3.10/site-packages/streamd/tests/ --run-md  
````
##### Test preparation functionality
````
pytest  Miniconda3/envs/md/lib/python3.10/site-packages/streamd/tests/ --run-preparation  
````

##### Test MD analysis functionality
````
pytest  Miniconda3/envs/md/lib/python3.10/site-packages/streamd/tests/ --run-analysis  
````
##### Test GBSA functionality
````
pytest  Miniconda3/envs/md/lib/python3.10/site-packages/streamd/tests/ --run-gbsa  
````
Test full functionality (will run tests for different cases): 
````
pytest  Miniconda3/envs/md/lib/python3.10/site-packages/streamd/tests/ --run-gbsa-full  
````
##### Test ProLIF functionality
````
pytest  Miniconda3/envs/md/lib/python3.10/site-packages/streamd/tests/ --run-prolif  
````
Test full functionality (will run tests for different cases): 
````
pytest  Miniconda3/envs/md/lib/python3.10/site-packages/streamd/tests/ --run-prolif-full  
````




