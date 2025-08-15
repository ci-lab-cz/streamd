### Clone repo
`git clone https://github.com/ci-lab-cz/streamd.git`

### Python Tests

To install pytest:
````
pip install -U pytest
````

#### Run tests:  
Add the ``--not-cleanup`` argument to ensure that test directories and test files are not removed and can be used for debugging purposes.
Add the ``-v`` or `-vs` arguments for printing verbose information.

##### Test full _md_run_ pipeline
````
pytest streamd/streamd/tests/ --run-md  
````
##### Test preparation functionality
````
pytest  streamd/streamd/tests/ --run-preparation  
````

##### Test MD analysis functionality
````
pytest streamd/streamd/tests/ --run-analysis  
````
##### Test GBSA functionality
````
pytest streamd/streamd/tests/ --run-gbsa  
````
Test full functionality (will run tests for different cases): 
````
pytest streamd/streamd/tests/ --run-gbsa-full  
````
##### Test ProLIF functionality
````
pytest streamd/streamd/tests/ --run-prolif  
````
Test full functionality (will run tests for different cases): 
````
pytest streamd/streamd/tests/ --run-prolif-full  
````

##### Test all
All the provided arguments above can be used together to run all tests in 1 pytest run:
````
pytest streamd/streamd/tests/  --run-preparation --run-analysis --run-gbsa-full --run-prolif-full --run-md
````




