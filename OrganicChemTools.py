# Import necessary libraries from RDKit for cheminformatics tasks
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
# Import matplotlib for image display
import matplotlib.pyplot as plt

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

def main():
    while True:
        # Main menu for the application
        print("\nOrganic Chemistry Assistant")
        print("1. Draw Chemical Structure")
        print("2. Calculate Molecular Weight")
        print("3. Exit")
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
            # Exit the program
            break
        else:
            # Handle invalid choices
            print("Invalid choice. Please try again.")

# Checks if the script is being run directly and not imported
if __name__ == "__main__":
    main()
