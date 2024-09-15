# Column-Design

The following Python script implements an iterative method using Numpy and Scipy libraries to determine the optimal column diamter to treat wastewater per client specifications.

The following equation is used for terminal velocity based on a force balance on a single particle in a fluidized bed:

<img width="150" alt="image" src="https://github.com/user-attachments/assets/64f12c34-38e1-4f14-aa5f-8f5a7c575f40">

Along with the following equation for the Reynolds number of a fluidized particle:

<img width="133" alt="image" src="https://github.com/user-attachments/assets/dbfa1ef8-48d1-4bf0-908a-bd6f315bafd1">

The drag coefficient was determined using the following literature correlations:

<img width="400" alt="image" src="https://github.com/user-attachments/assets/a483e6dc-47d5-45af-9256-e777ad1be3bf">

Once the terminal velocity is determined, the actual particle velocity can be optimized using an iterative process and the equations above to determine the optimal column diameter using the Richardson-Zaki equation:

<img width="100" alt="image" src="https://github.com/user-attachments/assets/d89bc4f1-b03e-406e-a9eb-d8416a35d675">

