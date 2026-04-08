

## Setup
This project uses Mojo for the simulation code and Python for visualization/analysis.

The recommended installation method is **Pixi**, which creates the full environment from [pixi.toml](pixi.toml), including:
- Mojo
- Python
- plotting/data libraries
- R dependencies used in notebooks

  Keep in mind, that currently (April 2026), mojo could be run by pixi only for Mac arm and Linux. To run it through other system, check Mojo installation guides

### Install with Pixi

1. Install Pixi
2. Clone this repository
3. Create/use the environment from `pixi.toml`

```bash
pixi install
pixi shell

### Running
To run each experiment enter the folder and run main.mojo file (recommended to look at it before running). The output of each experiment is csv file (visualisation is possible through uncommenting create_animated_visualization(simulation_data), but it was not used in the project).


Parameter Table: 

<img width="602" height="481" alt="parameter_table" src="https://github.com/user-attachments/assets/c50fdf65-b4ec-4e7c-bc59-1c4c50b41f29" />
