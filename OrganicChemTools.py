# Import necessary libraries from RDKit for cheminformatics tasks
import tkinter as tk
from tkinter import simpledialog
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, Crippen
# Import matplotlib for image display
import matplotlib.pyplot as plt
from PIL import Image, ImageTk

def draw_structure(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:  # Check if the molecule conversion was successful
            raise ValueError("Invalid SMILES notation.")
        Draw.MolToFile(mol, f"{smiles}.png")
        plt.imshow(plt.imread(f"{smiles}.png"))
        plt.axis('off')
        plt.show()
    except ValueError as e:
        print(e)
        return False  # Indicates failure
    return True  # Indicates success

def calculate_molecular_weight(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES notation.")
        return Descriptors.MolWt(mol)
    except ValueError as e:
        print(e)
        return None

# New function to calculate the logP value of a molecule
def calculate_logp(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:  # Check if the molecule conversion was successful
            raise ValueError("Invalid SMILES notation.")
        return Crippen.MolLogP(mol)
    except ValueError as e:
        print(e)
        return None

def on_draw_structure():
    smiles = simpledialog.askstring("Input", "Enter SMILES notation:", parent=root)
    if smiles:
        draw_structure(smiles)

def on_calculate_mw():
    smiles = simpledialog.askstring("Input", "Enter SMILES notation:", parent=root)
    if smiles:
        weight = calculate_molecular_weight(smiles)
        if weight is not None:
            result_label.config(text=f"Molecular Weight: {weight}")

def on_calculate_logp():
    smiles = simpledialog.askstring("Input", "Enter SMILES notation:", parent=root)
    if smiles:
        logp = calculate_logp(smiles)
        if logp is not None:
            result_label.config(text=f"LogP: {logp}")

# Create the main window
root = tk.Tk()
root.title("Organic Chemistry Assistant")

# Create buttons
draw_button = tk.Button(root, text="Draw Chemical Structure", command=on_draw_structure)
draw_button.pack()

mw_button = tk.Button(root, text="Calculate Molecular Weight", command=on_calculate_mw)
mw_button.pack()

logp_button = tk.Button(root, text="Calculate LogP", command=on_calculate_logp)
logp_button.pack()

# Label for displaying results
result_label = tk.Label(root, text="")
result_label.pack()

# Start the GUI loop
root.mainloop()

def main():
    while True:
        # Main menu for the application
        print("\nOrganic Chemistry Assistant")
        print("1. Draw Chemical Structure")
        print("2. Calculate Molecular Weight")
        print("3. Calculate LogP")
        print("4. Exit")
        choice = input("Enter your choice: ")

        # Handling the user's choice
        if choice == '1':
            while True:
                smiles = input("Enter SMILES notation: ")
                if draw_structure(smiles):
                    break  # Exit loop if successful
                else:
                    print("Please try again with a valid SMILES notation.")
        elif choice == '2':
            while True:
                smiles = input("Enter SMILES notation: ")
                weight = calculate_molecular_weight(smiles)
                if weight is not None:
                    print(f"Molecular Weight: {weight}")
                    break  # Exit loop if successful
                else:
                    print("Please try again with a valid SMILES notation.")
        elif choice == '3':
            while True:
                smiles = input("Enter SMILES notation: ")
                logp = calculate_logp(smiles)
                if logp is not None:
                    print(f"LogP: {logp}")
                    break  # Exit loop if successful
                else:
                    print("Please try again with a valid SMILES notation.")
        elif choice == '4':
            # Exit the program
            break
        else:
            # Handle invalid choices
            print("Invalid choice. Please try again.")

# Checks if the script is being run directly and not imported
if __name__ == "__main__":
    main()
