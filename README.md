## Setup
This project uses Mojo for the simulation code and Python for visualization and analysis.

The recommended installation method is Pixi, which creates the full environment from [pixi.toml](pixi.toml), including:
- Mojo
- Python
- plotting and data libraries
- R dependencies used in notebooks

Keep in mind that, as of April 2026, Mojo via Pixi is supported on macOS ARM and Linux. For other systems, check the official Mojo installation guide.

### Install with Pixi

1. Install Pixi
2. Clone this repository
3. Create and enter the environment

```bash
pixi install
pixi shell
```

### Running

To run an experiment, enter a particular experimtn folder and run the main Mojo file.

Each experiment writes CSV output. Visualization can also be enabled by uncommenting the visualization call in main.mojo. but this was not used for the main project results.
