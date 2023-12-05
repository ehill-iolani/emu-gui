# emu-gui
A shiny app for running the EMU pipeline through a GUI. This app also includes visualization tools for the EMU output.

## Installation
Build the docker image with the following commands:
```
git clone https://github.com/ehill-iolani/emu-gui.git
docker build -t emu-gui .
```

Pull the docker image from Docker Hub:
```
docker pull ethill/emu-gui:latest
```