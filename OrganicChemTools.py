# Import necessary libraries from RDKit for cheminformatics tasks
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
# Import matplotlib for image display
import matplotlib.pyplot as plt

def draw_structure(smiles):
    # Convert SMILES string to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    # Save the molecule's image as a PNG file using the SMILES string as the filename
    Draw.MolToFile(mol, f"{smiles}.png")
    # Display the image in a matplotlib plot
    plt.imshow(plt.imread(f"{smiles}.png"))
    plt.axis('off')  # Hide the axis for clarity
    plt.show()  # Show the plot with the molecule image

def calculate_molecular_weight(smiles):
    # Convert SMILES string to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    # Calculate and return the molecular weight of the molecule
    return Descriptors.MolWt(mol)

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
            # Draw chemical structure based on SMILES input
            smiles = input("Enter SMILES notation: ")
            draw_structure(smiles)
        elif choice == '2':
            # Calculate molecular weight based on SMILES input
            smiles = input("Enter SMILES notation: ")
            weight = calculate_molecular_weight(smiles)
            print(f"Molecular Weight: {weight}")
        elif choice == '3':
            # Exit the program
            break
        else:
            # Handle invalid choices
            print("Invalid choice. Please try again.")

# Checks if the script is being run directly and not imported
if __name__ == "__main__":
    main()
