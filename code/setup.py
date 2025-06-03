import subprocess
import os
import time
import venv

# Define the string to contain the splash screen welcome message
splash_screen_text: str = """ 
  _______________/ ____/ //_//   |     / /   /   |  / __ )_______________
 /____/____/____/ / __/ ,<  / /| |    / /   / /| | / __  /____/____/____/
/____/____/____/ /_/ / /| |/ ___ |   / /___/ ___ |/ /_/ /____/____/____/ 
    __    _____\____/_/ |_/_/__|_| _/_____/_/  |_/_____/_____________    
   / /   /  _/ ____/ / / /_  __/  / /   / __ \/ ____/ ____/ ____/ __ \   
  / /    / // / __/ /_/ / / /    / /   / / / / / __/ / __/ __/ / /_/ /   
 / /____/ // /_/ / __  / / /    / /___/ /_/ / /_/ / /_/ / /___/ _, _/    
/_____/___/\____/_/ /_/ /_/    /_____/\____/\____/\____/_____/_/ |_|                                                
"""


def main():
    # Display the startup message 
    for line in splash_screen_text.split('\n'):
        print(line)
        time.sleep(0.05)

    # Display for the user a list of potential environments that can be created
    requirements_dir: str =  os.path.join(os.path.dirname(__file__), "requirements")
    env_types: list[str] = os.listdir(requirements_dir)

    # Prompt the user to select a type 
    selected_type: None | str = None 
    while(selected_type is None):
        print("Please select the type of environment you would like to create:")
        for env_type_num, env_type in enumerate(env_types):
            print(f"\t{env_type_num} | {env_type.split('_')[0]}")

        # Retrieve the user's selection
        selection_index: str = input("").strip()

        # Validate the selection 
        if(not selection_index.isnumeric() or (int(selection_index) < 0 or int(selection_index) >= len(env_types))):
            print(f"ERROR: Invalid choice. Please choose a numeric valid between 0 and {len(env_types)}")
            continue
        
        # Otherwise, the selection was valid, and we can convert to int, get the type, and move forward 
        selected_type = env_types[int(selection_index)]

    # Retrieve the path to generate the virtual environment (will default to ./ if they just press enter)
    env_output_path: str = os.path.join(input("Please enter a path to save your environment in [./]: ").strip() or "./", f"{selected_type.split('_')[0]}_env")

    # Generate the environment
    print(f"Generating environment @ path: {env_output_path}")
    venv.create(env_output_path, with_pip=True)
    print("Generated!")
    
    # Install the relevant libraries
    print("Installing libraries...")

    # Retrieve the Python executable from the virtual environment 
    venv_python_executable: str = os.path.join(env_output_path, "bin", "python3")
    assert(os.path.exists(venv_python_executable))

    # Open the relevant requirements file
    with open(os.path.join(requirements_dir, selected_type), 'r') as f:
        # Iterate over the packages in the file 
        for line in f:
            # Retrieve the package from the line 
            package: str = line.strip()

            # If it was an empty line, we skip
            if(not package): continue

            # Otherwise, we will install into the desired virtual environment 
            subprocess.run([venv_python_executable, "-m", "pip", "install", package])

if(__name__ == "__main__"):
    main()