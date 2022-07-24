# Overview: 
A basic ORF(open reading frame) finder and DNA sequence analyzer.

## Getting started:
firstly clone the repository and ```cd``` into it and then later install all the required libraries with pip
```bash
pip3 install -r requirements.txt
```
and then later compile the shared object file
```bash
g++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) src/cpplib.cpp -o src/cpplib$(python3-config --extension-suffix)
```

then run the main.py in the src directory.
```bash
python3 src/main.py
```

### Technologies used
- Python 3.10
- C++11
- pybind11