import tkinter as tk
from tkinter import ttk, simpledialog, messagebox, scrolledtext
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, Crippen, AllChem, BRICS, rdMolDescriptors
from PIL import Image, ImageTk
import sqlite3
import logging
from typing import Optional

# --------------------------------------------------------------------------------
# Logging Setup
# --------------------------------------------------------------------------------
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# --------------------------------------------------------------------------------
# Application Configuration
# --------------------------------------------------------------------------------
class AppConfig:
    """
    A central place to define configuration parameters for the application.
    Using immutable class-level attributes to avoid mutable default pitfalls.
    """
    WINDOW_GEOMETRY = "900x700"
    DATABASE_NAME = "chemistry_data.db"
    CANVAS_SIZE = (400, 300)
    FUNCTIONAL_GROUPS = {
        "Alcohol": "[OH]",
        "Carboxylic Acid": "[CX3](=O)[OX2H1]",
        "Amine": "[NX3;H2,H1;!$(NC=O)]",
        "Carbonyl": "[CX3]=[OX1]",
        "Alkene": "[CX3]=[CX3]",
        "Alkyne": "[CX2]#[CX2]",
        "Ether": "[OX2]([CX4])[CX4]",
        "Ester": "[CX3](=O)[OX2][CX4]"
    }
    REACTION_TYPES = ["Oxidation", "Reduction", "Substitution", "Addition", "Elimination"]

# --------------------------------------------------------------------------------
# Utility Function: Generate 3D Structure
# --------------------------------------------------------------------------------
def generate_3d_structure(mol: Optional[Chem.Mol]) -> Optional[Chem.Mol]:
    """
    Generate a 3D conformation for an RDKit molecule.
    Returns a new molecule with 3D coordinates embedded and optimized.

    Parameters:
        mol (rdkit.Chem.Mol): The input RDKit Mol object (2D or partial).
    
    Returns:
        rdkit.Chem.Mol: A new Mol object with 3D coordinates or None if generation fails.
    """
    if mol is None:
        return None
    try:
        # Add hydrogens to improve 3D embedding
        mol_3d = Chem.AddHs(mol)
        # Embed 3D coordinates
        AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())
        # Optimize geometry
        AllChem.MMFFOptimizeMolecule(mol_3d)
        return mol_3d
    except Exception as e:
        logger.error(f"3D generation failed: {e}")
        return None

# --------------------------------------------------------------------------------
# Main Application Class
# --------------------------------------------------------------------------------
class OrganicChemistryAssistant:
    """
    A Tkinter-based GUI for advanced organic chemistry assistance.
    """
    def __init__(self, root: tk.Tk):
        """
        Initialize the GUI components, create the database, and set up the tabs.
        """
        self.root = root
        self.root.title("Advanced Organic Chemistry Assistant")
        self.root.geometry(AppConfig.WINDOW_GEOMETRY)

        # Database setup
        self.conn = sqlite3.connect(AppConfig.DATABASE_NAME)
        self.cursor = self.conn.cursor()
        self._create_database_table()

        # Create notebook for tabs
        self.notebook = ttk.Notebook(root)
        self.notebook.pack(expand=True, fill='both', padx=5, pady=5)

        # Create the tabs
        self.structure_tab = self._create_tab("Structure")
        self.properties_tab = self._create_tab("Properties")
        self.analysis_tab = self._create_tab("Analysis")
        self.reaction_tab = self._create_tab("Reactions")
        self.database_tab = self._create_tab("Database")

        # Setup tab contents
        self.setup_structure_tab()
        self.setup_properties_tab()
        self.setup_analysis_tab()
        self.setup_reaction_tab()
        self.setup_database_tab()

        # Status bar
        self.status_var = tk.StringVar()
        self.status_bar = ttk.Label(root, textvariable=self.status_var, relief=tk.SUNKEN)
        self.status_bar.pack(side=tk.BOTTOM, fill=tk.X)

        # Internal image reference (to prevent garbage collection in Tk)
        self.photo = None

    def _create_database_table(self) -> None:
        """
        Create the necessary database table for storing molecules if it does not exist.
        """
        self.cursor.execute("""
            CREATE TABLE IF NOT EXISTS Molecules (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                smiles TEXT NOT NULL,
                properties TEXT
            )
        """)

    def _create_tab(self, title: str) -> ttk.Frame:
        """
        Helper function to create a tab within the notebook.
        """
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text=title)
        return tab

    # ----------------------------------------------------------------------------
    # Structure Tab
    # ----------------------------------------------------------------------------
    def setup_structure_tab(self):
        """
        Create the layout and widgets for the 'Structure' tab, including the SMILES entry,
        draw/clear buttons, and a canvas to display the molecule image.
        """
        frame = ttk.LabelFrame(self.structure_tab, text="Molecule Visualization")
        frame.pack(expand=True, fill='both', padx=5, pady=5)

        # SMILES input
        ttk.Label(frame, text="SMILES:").pack(padx=5, pady=2)
        self.smiles_entry = ttk.Entry(frame, width=50)
        self.smiles_entry.pack(padx=5, pady=2)

        # Action buttons
        btn_frame = ttk.Frame(frame)
        btn_frame.pack(pady=5)

        ttk.Button(btn_frame, text="Draw Structure", command=self.draw_structure).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Clear", command=self.clear_structure).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Generate 3D", command=self.on_generate_3d).pack(side=tk.LEFT, padx=5)

        # Canvas for molecule display
        self.canvas = tk.Canvas(frame, bg='white',
                                width=AppConfig.CANVAS_SIZE[0],
                                height=AppConfig.CANVAS_SIZE[1])
        self.canvas.pack(expand=True, fill='both', padx=5, pady=5)

    def draw_structure(self):
        """
        Draw the 2D structure of the molecule from the SMILES string on the canvas.
        """
        smiles = self.smiles_entry.get().strip()
        if not smiles:
            messagebox.showerror("Error", "Please enter a SMILES string.")
            return
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError("Invalid SMILES notation")

            img = Draw.MolToImage(mol)
            self.photo = ImageTk.PhotoImage(img)
            self.canvas.delete("all")
            self.canvas.create_image(AppConfig.CANVAS_SIZE[0] // 2,
                                     AppConfig.CANVAS_SIZE[1] // 2,
                                     image=self.photo)
            self.status_var.set("Structure drawn successfully.")
        except Exception as e:
            logger.error("Error drawing structure: %s", e)
            messagebox.showerror("Error", str(e))
            self.status_var.set("Error drawing structure.")

    def clear_structure(self):
        """
        Clear the canvas and the SMILES entry box.
        """
        self.canvas.delete("all")
        self.smiles_entry.delete(0, tk.END)
        self.status_var.set("Canvas cleared.")

    def on_generate_3d(self):
        """
        Generate a 3D structure from the SMILES input, then display a 2D depiction
        of the 3D-optimized molecule on the canvas.
        """
        smiles = self.smiles_entry.get().strip()
        if not smiles:
            messagebox.showerror("Error", "Please enter a SMILES string.")
            return

        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                raise ValueError("Invalid SMILES notation")
            mol_3d = generate_3d_structure(mol)

            if mol_3d is None:
                raise ValueError("3D generation returned None.")

            # Display 2D depiction of the 3D-optimized molecule
            img = Draw.MolToImage(mol_3d)
            self.photo = ImageTk.PhotoImage(img)
            self.canvas.delete("all")
            self.canvas.create_image(AppConfig.CANVAS_SIZE[0] // 2,
                                     AppConfig.CANVAS_SIZE[1] // 2,
                                     image=self.photo)
            self.status_var.set("3D structure generated and displayed (2D depiction).")
        except Exception as e:
            logger.error("Error generating 3D structure: %s", e)
            messagebox.showerror("Error", str(e))
            self.status_var.set("Error generating 3D structure.")

    # ----------------------------------------------------------------------------
    # Properties Tab
    # ----------------------------------------------------------------------------
    def setup_properties_tab(self):
        """
        Create the layout for the 'Properties' tab, including a textbox to display
        calculated properties and a button to trigger property calculations.
        """
        frame = ttk.LabelFrame(self.properties_tab, text="Molecular Properties")
        frame.pack(expand=True, fill='both', padx=5, pady=5)

        self.properties_text = scrolledtext.ScrolledText(frame, height=20, width=60)
        self.properties_text.pack(padx=5, pady=5)

        ttk.Button(frame, text="Calculate Properties", command=self.calculate_properties).pack(pady=5)

    def calculate_properties(self):
        """
        Calculate and display various molecular properties based on the SMILES input,
        including molecular weight, logP, TPSA, number of rotatable bonds, etc.
        """
        smiles = self.smiles_entry.get().strip()
        if not smiles:
            messagebox.showerror("Error", "Please enter a SMILES string.")
            return
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError("Invalid SMILES notation")

            properties = {
                "Molecular Formula": rdMolDescriptors.CalcMolFormula(mol),
                "Molecular Weight": Descriptors.ExactMolWt(mol),
                "LogP (Crippen)": Crippen.MolLogP(mol),
                "TPSA": Descriptors.TPSA(mol),
                "Number of Rotatable Bonds": Descriptors.NumRotatableBonds(mol),
                "Number of H-Bond Donors": Descriptors.NumHDonors(mol),
                "Number of H-Bond Acceptors": Descriptors.NumHAcceptors(mol),
                "Number of Rings": rdMolDescriptors.CalcNumRings(mol),
            }

            self.properties_text.delete(1.0, tk.END)
            for prop, value in properties.items():
                # If the value is float-like, format it with two decimals
                if isinstance(value, float):
                    self.properties_text.insert(tk.END, f"{prop}: {value:.2f}\n")
                else:
                    self.properties_text.insert(tk.END, f"{prop}: {value}\n")

            self.status_var.set("Properties calculated successfully.")
        except Exception as e:
            logger.error("Error calculating properties: %s", e)
            messagebox.showerror("Error", str(e))
            self.status_var.set("Error calculating properties.")

    # ----------------------------------------------------------------------------
    # Analysis Tab
    # ----------------------------------------------------------------------------
    def setup_analysis_tab(self):
        """
        Create the layout for the 'Analysis' tab, including radio buttons to select
        fragment or functional group analysis, and a textbox for displaying results.
        """
        frame = ttk.LabelFrame(self.analysis_tab, text="Structure Analysis")
        frame.pack(expand=True, fill='both', padx=5, pady=5)

        self.analysis_var = tk.StringVar(value="fragments")
        ttk.Radiobutton(frame, text="Fragment Analysis", variable=self.analysis_var, value="fragments").pack()
        ttk.Radiobutton(frame, text="Functional Groups", variable=self.analysis_var, value="groups").pack()

        ttk.Button(frame, text="Analyze", command=self.analyze_structure).pack(pady=5)

        self.analysis_text = scrolledtext.ScrolledText(frame, height=15, width=60)
        self.analysis_text.pack(padx=5, pady=5)

    def analyze_structure(self):
        """
        Perform structure analysis based on the selected option (fragments or
        functional groups), then display the results in the analysis text box.
        """
        smiles = self.smiles_entry.get().strip()
        if not smiles:
            messagebox.showerror("Error", "Please enter a SMILES string.")
            return
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError("Invalid SMILES notation")

            analysis_type = self.analysis_var.get()
            self.analysis_text.delete(1.0, tk.END)

            if analysis_type == "fragments":
                fragments = list(BRICS.BRICSDecompose(mol))
                self.analysis_text.insert(tk.END, "BRICS Fragments:\n\n")
                for i, frag in enumerate(fragments, 1):
                    self.analysis_text.insert(tk.END, f"{i}. {frag}\n")
            elif analysis_type == "groups":
                self.analysis_text.insert(tk.END, "Functional Groups Found:\n\n")
                for group, smarts in AppConfig.FUNCTIONAL_GROUPS.items():
                    pattern = Chem.MolFromSmarts(smarts)
                    if mol.HasSubstructMatch(pattern):
                        matches = mol.GetSubstructMatches(pattern)
                        self.analysis_text.insert(tk.END, f"{group}: {len(matches)} occurrence(s)\n")

            self.status_var.set("Analysis completed successfully.")
        except Exception as e:
            logger.error("Error analyzing structure: %s", e)
            messagebox.showerror("Error", str(e))
            self.status_var.set("Error during analysis.")

    # ----------------------------------------------------------------------------
    # Reaction Tab
    # ----------------------------------------------------------------------------
    def setup_reaction_tab(self):
        """
        Create the layout for the 'Reactions' tab, including a dropdown to select
        the reaction type and a button to predict products.
        """
        frame = ttk.LabelFrame(self.reaction_tab, text="Reaction Prediction")
        frame.pack(expand=True, fill='both', padx=5, pady=5)

        ttk.Label(frame, text="Reaction Type:").pack()
        self.reaction_type = ttk.Combobox(frame, values=AppConfig.REACTION_TYPES)
        self.reaction_type.pack(pady=5)

        ttk.Button(frame, text="Predict Products", command=self.predict_reaction).pack(pady=5)

        self.reaction_text = scrolledtext.ScrolledText(frame, height=15, width=60)
        self.reaction_text.pack(padx=5, pady=5)

    def predict_reaction(self):
        """
        Predict possible reaction outcomes based on the chosen reaction type
        and the current SMILES string.
        """
        smiles = self.smiles_entry.get().strip()
        if not smiles:
            messagebox.showerror("Error", "Please enter a SMILES string.")
            return
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError("Invalid SMILES notation")

            reaction_type = self.reaction_type.get()
            self.reaction_text.delete(1.0, tk.END)

            if reaction_type == "Oxidation":
                if mol.HasSubstructMatch(Chem.MolFromSmarts("[CH2][OH]")):
                    self.reaction_text.insert(tk.END, "Predicted: Alcohol oxidation to aldehyde/ketone\n")
                elif mol.HasSubstructMatch(Chem.MolFromSmarts("[CH3]")):
                    self.reaction_text.insert(tk.END, "Predicted: Alkane oxidation possible\n")
            elif reaction_type == "Reduction":
                if mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3]=[OX1]")):
                    self.reaction_text.insert(tk.END, "Predicted: Carbonyl reduction to alcohol\n")
                elif mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3]=[CX3]")):
                    self.reaction_text.insert(tk.END, "Predicted: Alkene reduction to alkane\n")
            elif reaction_type == "Substitution":
                if mol.HasSubstructMatch(Chem.MolFromSmarts("[CX4][Br,I,Cl,F]")):
                    self.reaction_text.insert(tk.END, "Predicted: Halide substitution possible\n")
            elif reaction_type == "Addition":
                if mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3]=[CX3]")):
                    self.reaction_text.insert(tk.END, "Predicted: Alkene addition possible\n")
            elif reaction_type == "Elimination":
                if mol.HasSubstructMatch(Chem.MolFromSmarts("[CX4][OH]")):
                    self.reaction_text.insert(tk.END, "Predicted: Alcohol elimination to alkene\n")
            else:
                self.reaction_text.insert(tk.END, "Reaction type not recognized or not implemented\n")

            self.status_var.set("Reaction prediction completed.")
        except Exception as e:
            logger.error("Error predicting reaction: %s", e)
            messagebox.showerror("Error", str(e))
            self.status_var.set("Error predicting reaction.")

    # ----------------------------------------------------------------------------
    # Database Tab
    # ----------------------------------------------------------------------------
    def setup_database_tab(self):
        """
        Create the layout for the 'Database' tab, including buttons to save
        the current molecule and retrieve all saved data.
        """
        frame = ttk.LabelFrame(self.database_tab, text="Database Management")
        frame.pack(expand=True, fill='both', padx=5, pady=5)

        ttk.Button(frame, text="Save Molecule", command=self.save_molecule).pack(pady=5)
        ttk.Button(frame, text="Retrieve Data", command=self.retrieve_data).pack(pady=5)

        self.database_text = scrolledtext.ScrolledText(frame, height=20, width=80)
        self.database_text.pack(padx=5, pady=5)

    def save_molecule(self):
        """
        Save the current SMILES and property data to the database table.
        """
        smiles = self.smiles_entry.get().strip()
        properties = self.properties_text.get(1.0, tk.END).strip()

        if not smiles or not properties:
            messagebox.showerror("Error", "SMILES or properties are empty!")
            return

        try:
            self.cursor.execute("INSERT INTO Molecules (smiles, properties) VALUES (?, ?)",
                                (smiles, properties))
            self.conn.commit()
            messagebox.showinfo("Success", "Molecule saved successfully.")
            self.status_var.set("Molecule saved to database.")
        except Exception as e:
            logger.error("Error saving molecule: %s", e)
            messagebox.showerror("Error", "Error saving molecule: " + str(e))

    def retrieve_data(self):
        """
        Retrieve all molecule records from the database and display them
        in the text box.
        """
        try:
            self.cursor.execute("SELECT * FROM Molecules")
            rows = self.cursor.fetchall()

            self.database_text.delete(1.0, tk.END)
            for row in rows:
                self.database_text.insert(tk.END, f"ID: {row[0]}\nSMILES: {row[1]}\nProperties:\n{row[2]}\n")
                self.database_text.insert(tk.END, f"{'-'*50}\n")
            self.status_var.set("Data retrieved from database.")
        except Exception as e:
            logger.error("Error retrieving data: %s", e)
            messagebox.showerror("Error", "Error retrieving data: " + str(e))

# --------------------------------------------------------------------------------
# Main Execution
# --------------------------------------------------------------------------------
if __name__ == "__main__":
    root = tk.Tk()
    app = OrganicChemistryAssistant(root)
    root.mainloop()
