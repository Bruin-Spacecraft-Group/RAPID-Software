# Rapid-0-Software

**Rapid-0-Software** contains all of the flight and testing software for RAPID-0. 

Refer to the [How-to Software](https://docs.google.com/document/d/1KSEFYbrn3GWrpJz5_CZqGYSEJwLImQULD9oqK3tNjNo/edit?tab=t.0) document for more information related to the project RAPID Software Subteam practices


## Repository Structure:

- .github: Contains information about running tests on push and owners of the repo
- artifacts: The software deployed to a specific board
    - board.py contains all the pin definitions used by that artifact (except if it’s deployed to the Nucleo board, that board uses its own board.py file that overwrites the one written here)
    - code.py is the main part of the artifact, it is the code that will be running on the artifact during operation
    - include.json contains a list of the dependencies required for the artifact from the src folder and/or the unit_tests folder
- config: Sets up the running of unit tests
- src: Contains the bulk of the code, including the tasks and code to run the tasks
    - drivers contains the code necessary for interfacing with hardware elements
    - lib contains libraries that don’t directly interact with the hardware, and instead exist only as software
    - tasks contains each individual application with a specific purpose, and are the tasks that will be performed in an artifact. Tasks can be used for one artifact only, or could be applied to multiple different artifacts
- submodules: contains modules that rely on circuitpython but aren’t included in base circuitpython
- tools: contains scripts useful for code development and code application
    - deploy_to_usb is how to transfer software onto a board via usb cable
- unit_tests: tests to check that a library or task is working properly

## Cloning the Repository

In order to create a local copy of the RAPID-0 Software repository, `clone` the repo (which will create a copy in your current directory):

    git clone https://github.com/Bruin-Spacecraft-Group/RAPID-0-Software.git
        
After cloning the repository, be sure to run

    git submodule update --init <path to submodule>

for each submodule (`path to submodule` usually will just be `submodules/Adafruit_CircuitPython_asyncio` and `submodules/Adafruit_CircuitPython_Ticks`), this updates the submodules (other git repos that have dependencies we use)
