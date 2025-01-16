# IMPORTS
import sys
import os
import pandas as pd
import re
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPixmap, QFont, QIcon, QClipboard, QBrush, QColor
from PyQt5.QtWidgets import QApplication, QMainWindow, QListWidget, QLineEdit, QPushButton, QMessageBox, QLabel, QFileDialog, QRadioButton
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# AUTHOR - Shean Mobed
# ORG - Biosurv International

# DESCRIPTION
"""

"""

# WIN Compile Command
"""
nuitka --onefile --enable-plugins=pyqt5 --include-data-files=Logo.png=./Logo.png --include-data-files=Icon.ico=./Icon.ico --disable-console 
--windows-icon-from-ico=Icon.ico --company-name="Biosurv International" 
--product-name="FASTAQC" --file-version=1.0.0 
--file-description=="This App curates Fasta seqs to remove sequences that failed DDNS QC"  
"""

# BUILD PROGRESS
"""
- make window with directory selection for piranha report 
- make text box or selection box for sequences, perhaps tabs for sabins for quicker selection 
- build logic that reads report and pre selects failed sequences 
- save button will save to vp1 sequences with _QCPass tag 
- make double click feature for failing samples in app, disabled in preference for in report changes - disabled feature, code commented out line 434 - 464
- add option to choose between ddns and isolat 
- add checks for is samples match between fasta and report, empty QC cells, if they try to generate fasta before parsing, file presence and handling checks 

- show drop box bg label when done with testing

- added version number in bottom righthand corner 
- now removes non-polio seqs from output and non samples based of name regex
- generate buttton is disabled until files parsed, colour set to darkgray
- checks RunQC to fail samples
"""

# APP
class ListBoxWidget(QListWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAcceptDrops(True)
        self.resize(500, 500)


    def dragEnterEvent(self, event):
        #bg_label.hide()
        if event.mimeData().hasUrls():
            urls = event.mimeData().urls()
            if all(url.toString().endswith(('.csv','.fasta')) for url in urls) and len(urls) <= 2:
                event.accept()
            else:
                event.ignore()
        else:
            event.ignore()

    def dragMoveEvent(self, event):
        if event.mimeData().hasUrls():
            event.setDropAction(Qt.CopyAction)
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event):
        #self.clear()
        urls = event.mimeData().urls()
        file_paths = [url.toLocalFile() for url in urls if url.toString().endswith(('.csv','.fasta'))]
        if self.count() < 2:
            self.addItems(file_paths)
        else:
            QMessageBox.warning(self, 'Warning', 'No more than two files can be provided')
            return
        #bg_label.setText('')

class SelectionListBox(QListWidget):
    pass
    

class App(QMainWindow):
    def __init__(self):
        super().__init__()

        screen_size = self.screen().size()
        screen_width = int(screen_size.width() * 0.5)
        screen_height = int(screen_size.height() * 0.8)
        self.resize(screen_width, screen_height)  # 1200,1000 original

        self.setStyleSheet("background-color: white; color: black;")  # Set window color to white

        self.setWindowIcon(
            QIcon(os.path.join(os.path.dirname(__file__), "psc_logo.ico")))  # Sets top left icon to custom icon
        self.setWindowTitle("Fasta QC")  # App Title
        
        # ---- VERSION ---- # change when new veriosn
        self.app_ver = QLabel('Version: 1.0.0',self)
        self.app_ver.setStyleSheet("background-color:transparent")
        self.app_ver.setGeometry(int(screen_width * 0.88), int(screen_height * 0.96), int(screen_width * 0.2), int(screen_height * 0.06))
        self.app_ver.setFont(QFont('Arial', 9))
        
        # ---- LOGO ----
        self.logo_label = QLabel(self)
        self.logo_label.setGeometry(int(screen_width * 0.02), int(screen_height * 0.02), int(screen_width * 0.1), int(screen_height * 0.1))  # Logo position x y and dimension x y
        pixmap = QPixmap(os.path.join(os.path.dirname(__file__), 'psc_logo.png'))
        self.logo_label.setPixmap(pixmap)
        self.logo_label.setScaledContents(True)
        
        # ---- TITLE ----
        self.app_title = QLabel('Quality Controlled and Database Fasta Generator', self)
        self.app_title.setGeometry(int(screen_width * 0.2), int(screen_height * 0.02), int(screen_width * 0.7), int(screen_height * 0.1))  # Logo position x y and dimension x y
        self.app_title.setFont(QFont('Arial', 15))
        
        # ---- DESTINATION BOX ----
        self.destination_entry = QLineEdit(self)
        self.destination_entry.setText("C:/Users/SheanMobed/OneDrive - Biosurv International/Desktop") # TEMP
        self.destination_entry.setStyleSheet("background-color:#FAF9F6;border-color:lightgrey;border-style:dashed;border-width: 2px;border-radius: 10px;")
        self.destination_entry.setGeometry(int(screen_width * 0.13), int(screen_height * 0.27), int(screen_width * 0.8), int(screen_height * 0.04))
        self.destination_entry.setFont(QFont('Arial', 11))

        self.destination_label = QLabel('Destination:', self)
        self.destination_label.setStyleSheet("background-color:transparent")
        self.destination_label.setGeometry(int(screen_width * 0.03), int(screen_height * 0.265), int(screen_width * 0.2), int(screen_height * 0.05))
        self.destination_label.setFont(QFont('Arial', 9))

        # ---- DROPBOX ----
        self.listbox_view = ListBoxWidget(self)
        self.listbox_view.setStyleSheet("background-color:#FAF9F6;border-color:lightgrey;border-style:dashed;border-width: 2px;border-radius: 10px;")
        self.listbox_view.setGeometry(int(screen_width * 0.13), int(screen_height * 0.19), int(screen_width * 0.8), int(screen_height * 0.08))
        
        self.listbox_view.addItem('C:/Users/SheanMobed/Documents/Coding/Polio/DDNS_Polio/Final_Data/20240702_035_detailed_run_report.csv') # TEMP
        self.listbox_view.addItem('C:/Users/SheanMobed/Documents/Coding/Polio/DDNS_Polio/Miscellaneous/2024/20240702_035/20240702_035_vp1_sequences.fasta')
        
        self.listbox_label = QLabel('Dropbox:', self)
        self.listbox_label.setStyleSheet("background-color:transparent")
        self.listbox_label.setGeometry(int(screen_width * 0.05), int(screen_height * 0.21), int(screen_width * 0.2), int(screen_height * 0.02))
        self.listbox_label.setFont(QFont('Arial', 9))

        # global bg_label # set to global scope so it can be changed in various other app functions
        # bg_label = QLabel('Drop CSVs here', self)
        # bg_label.setStyleSheet("background-color:#FAF9F6")
        # bg_label.setGeometry(int(screen_width * 0.5), int(screen_height * 0.18), int(screen_width * 0.2), int(screen_height * 0.05))
        
        #---- S1 LIST ----
        self.s1_list = SelectionListBox(self)
        self.s1_list.setStyleSheet("background-color:#FAF9F6;border-color:lightgrey;border-width: 2px;border-radius: 10px;")
        self.s1_list.setGeometry(int(screen_width * 0.03), int(screen_height * 0.34), int(screen_width * 0.3), int(screen_height * 0.55))
        self.s1_list.setFont(QFont('Arial', 9))
        
        self.s1_label = QLabel('POLIO 1',self)
        self.s1_label.setStyleSheet("background-color:#2e3192; color: white; border-radius: 10px")
        self.s1_label.setGeometry(int(screen_width * 0.03), int(screen_height * 0.32), int(screen_width * 0.3), int(screen_height * 0.02))
        self.s1_label.setFont(QFont('Arial', 9))
        self.s1_label.setAlignment(Qt.AlignHCenter)
        
        
        #---- S2 LIST ----
        self.s2_list = SelectionListBox(self)
        self.s2_list.setStyleSheet("background-color:#FAF9F6;border-color:lightgrey;border-width: 2px;border-radius: 10px;")
        self.s2_list.setGeometry(int(screen_width * 0.35), int(screen_height * 0.34), int(screen_width * 0.3), int(screen_height * 0.55))
        self.s2_list.setFont(QFont('Arial', 9))
        
        self.s2_label = QLabel('POLIO 2',self)
        self.s2_label.setStyleSheet("background-color:#2e3192; color: white; border-radius: 10px")
        self.s2_label.setGeometry(int(screen_width * 0.35), int(screen_height * 0.32), int(screen_width * 0.3), int(screen_height * 0.02))
        self.s2_label.setFont(QFont('Arial', 9))
        self.s2_label.setAlignment(Qt.AlignHCenter)
        
        #---- S3 LIST ----
        self.s3_list = SelectionListBox(self)
        self.s3_list.setStyleSheet("background-color:#FAF9F6;border-color:lightgrey;border-width: 2px;border-radius: 10px;")
        self.s3_list.setGeometry(int(screen_width * 0.67), int(screen_height * 0.34), int(screen_width * 0.3),int(screen_height * 0.55))
        self.s3_list.setFont(QFont('Arial', 9))
        
        self.s3_label = QLabel('POLIO 3',self)
        self.s3_label.setStyleSheet("background-color:#2e3192; color: white; border-radius: 10px")
        self.s3_label.setGeometry(int(screen_width * 0.67), int(screen_height * 0.32), int(screen_width * 0.3), int(screen_height * 0.02))
        self.s3_label.setFont(QFont('Arial', 9))
        self.s3_label.setAlignment(Qt.AlignHCenter)
        
        
        # ---- BUTTONS ----
        # save button        
        self.btn_save = QPushButton('Generate \nQC Fasta', self)
        self.btn_save.setGeometry(int(screen_width * 0.54), int(screen_height * 0.9), int(screen_width * 0.15), int(screen_height * 0.07))  
        self.btn_save.setFont(QFont('Arial', 11))
        self.btn_save.setStyleSheet("QPushButton"
                                     "{"
                                     "color: white; border-radius: 15px;background-color:#2e3192;"
                                     "border-color:black;border-style: solid;border-width: 1px;"
                                     "}"
                                     "QPushButton::pressed"
                                     "{"
                                     "background-color : #3638d8;"
                                     "}"
                                     "QPushButton:disabled { background-color: grey; color: darkgrey; }")
        self.btn_save.setEnabled(False)
        self.btn_save.clicked.connect(self.save_fasta)
        
        # parse button
        self.btn_parse = QPushButton('Parse Files', self)
        self.btn_parse.setGeometry(int(screen_width * 0.31), int(screen_height * 0.9), int(screen_width * 0.15), int(screen_height * 0.07))  
        self.btn_parse.setFont(QFont('Arial', 11))
        self.btn_parse.setStyleSheet("QPushButton"
                                     "{"
                                     "color: white; border-radius: 15px;background-color:#2e3192;"
                                     "border-color:black;border-style: solid;border-width: 1px;"
                                     "}"
                                     "QPushButton::pressed"
                                     "{"
                                     "background-color : #3638d8;"
                                     "}")
        self.btn_parse.clicked.connect(self.parse_data)

        
        # User defined file destination
        self.destination_btn = QPushButton('Destination', self)
        self.destination_btn.setGeometry(int(screen_width * 0.08), int(screen_height * 0.9), int(screen_width * 0.15), int(screen_height * 0.07))  # 750, 770, 300, 100
        self.destination_btn.setFont(QFont('Arial', 11))
        self.destination_btn.setStyleSheet("QPushButton"
                                           "{"
                                           "color: white; border-radius: 15px;background-color:#2e3192;"
                                           "border-color:black;border-style: solid;border-width: 1px;"
                                           "}"
                                           "QPushButton::pressed"
                                           "{"
                                           "background-color : #3638d8;"
                                           "}")
        
        self.destination_btn.clicked.connect(lambda: self.select_destination(3))
        
        self.btn_clear = QPushButton('Clear', self)
        self.btn_clear.setGeometry(int(screen_width * 0.76), int(screen_height * 0.9), int(screen_width * 0.15), int(screen_height * 0.07))  
        self.btn_clear.setFont(QFont('Arial', 11))
        self.btn_clear.setStyleSheet("QPushButton"
                                     "{"
                                     "color: white; border-radius: 15px;background-color:#2e3192;"
                                     "border-color:black;border-style: solid;border-width: 1px;"
                                     "}"
                                     "QPushButton::pressed"
                                     "{"
                                     "background-color : #3638d8;"
                                     "}")
        self.btn_clear.clicked.connect(self.clear_list)
        
        
        # ---- METHOD RADIO BUTTONS ----
        self.radio_ddns = QRadioButton('DDNS', self)
        self.radio_ddns.setGeometry(int(screen_width * 0.14), int(screen_height * 0.15), int(screen_width * 0.1), int(screen_height * 0.02))
        self.radio_ddns.setChecked(True)

        self.radio_isolat = QRadioButton('Isolate', self)
        self.radio_isolat.setGeometry(int(screen_width * 0.22), int(screen_height * 0.15), int(screen_width * 0.1), int(screen_height * 0.02))

        self.radio_label = QLabel('Run Method:',self)
        self.radio_label.setGeometry(int(screen_width * 0.03), int(screen_height * 0.15), int(screen_width * 0.1), int(screen_height * 0.02))
        self.radio_label.setFont(QFont('Arial', 9))

        
    def select_destination(self, type):
        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.AnyFile if type in (1, 2) else QFileDialog.Directory)
        
        if file_dialog.exec_():
            selected_files = file_dialog.selectedFiles()
            if selected_files:
                if type == 1:
                    self.epi_entry.setText(selected_files[0])
                    self.listbox_view.addItem(selected_files[0])
                    #bg_label.hide()

                elif type == 2:
                    self.lab_entry.setText(selected_files[0])
                    self.listbox_view.addItem(selected_files[0])
                    #bg_label.hide()

                elif type == 3:
                    self.destination_entry.setText(selected_files[0])
                

    def parse_data(self):
        # handle inputs
        paths = set(self.listbox_view.item(i).text() for i in range(self.listbox_view.count()))
        paths = sorted(paths)
        
        # clear list to initialise
        self.s1_list.clear()
        self.s2_list.clear()
        self.s3_list.clear()
        
        if len(paths) == 0:
            QMessageBox.warning(self, 'Warning', 'No Files selected')
            return
        
        global destination_path
        destination_path = self.destination_entry.text()
        
        if destination_path == '':
            QMessageBox.warning(self, 'Warning', 'No destination selected')
            return
        try:
            report = [file for file  in paths if file.endswith('.csv')][0]
        except IndexError:
            QMessageBox.warning(self, 'Warning', 'No Report Given')
            return
        
        try:
            global fasta
            fasta  = [file for file  in paths if file.endswith('.fasta')][0]
        except IndexError:
            QMessageBox.warning(self, 'Warning', 'No Fasta Given')
            return
        
        # readin report
        report = pd.read_csv(report, sep=None)
        
        # Check if Fasta samples match report samples
        discordant_fasta_samples = []
        for sequence in SeqIO.parse(fasta, 'fasta'):
            
            # Extracts sample ID, corrects if '-DDNS' is present
            sample_id = sequence.id.split('|')[0].replace('-DDNS', '')
            # Check if sample is in the report
            if sample_id not in report['sample'].values:
                discordant_fasta_samples.append(sample_id)

        # warning if discordant samples are found
        if discordant_fasta_samples:
            sample_list = ', '.join(discordant_fasta_samples)
            QMessageBox.warning(self, 'Warning', f'The following samples are present in FASTA but not in the report: {sample_list}.\n Please check if they are from the same run.')
            return
        
        # Filter for QC and Sample Passes, fill empty cells with empty string to not cause attribute error
        report['RunQC'] = report['RunQC'].fillna('').str.strip().str.title()
        report['SampleQC'] = report['SampleQC'].fillna('').str.strip().str.title()
        
        # Check for empty strings in 'RunQC' and 'SampleQC' columns
        empty_qc_rows = report[(report['RunQC'] == '') | (report['SampleQC'] == '')]

        # If there are any empty strings, show a QMessageBox
        if not empty_qc_rows.empty:
            empty_samples = empty_qc_rows['sample'].tolist()  # Assuming 'sample' column identifies each sample
            sample_list = ', '.join(empty_samples)
            
            QMessageBox.warning(self, 'Warning', f'The following samples have empty values in RunQC or SampleQC:\n{sample_list}\nPlease check the report.')
            return

        # Standardising DDNS Classification for report
        def classify(row):
            classifications = []
            if row['Sabin1-related|classification'] == 'Sabin-like':
                classifications.append('SABIN1')
            elif row['Sabin1-related|classification'] == 'VDPV':
                classifications.append('VDPV1')

            if row['Sabin2-related|classification'] == 'Sabin-like':
                classifications.append('SABIN2')
            elif row['Sabin2-related|classification'] == 'VDPV':
                classifications.append('VDPV2')

            if row['Sabin3-related|classification'] == 'Sabin-like':
                classifications.append('SABIN3')
            elif row['Sabin3-related|classification'] == 'VDPV':
                classifications.append('VDPV3')
            
            if row['NonPolioEV|classification'] == 'NonPolioEV':
                classifications.append('NONPOLIOEV')

            return '+'.join(classifications) if classifications else 'Negative'

        # Apply the classify function to each row
        report['DDNSclassification'] = report.apply(classify, axis=1)
        
        # Remove controls from df using regex
        report = report.loc[report['sample'].str.match(r'^[A-Z]{3}-\d{2}-\d{3,5}$')]
        
        # Splitting on the + and create duplicate row below with explode
        if not report[report['DDNSclassification'].str.contains('\\+', na=False)].empty:
            combos = report[report['DDNSclassification'].str.contains('\\+', na=False)]
            combos.loc[:,'DDNSclassification'] = combos['DDNSclassification'].str.split("\\+")
            combos = combos.explode('DDNSclassification')
            report = pd.concat([report,combos])
        
        # Removes non polio after split
        report = report.loc[(report['DDNSclassification'] != 'NONPOLIOEV')]
        
        # sorted number for easier reading
        report = report.sort_values(by='sample')

        # PASS FAIL
        ddns_pass = report.loc[(report['SampleQC'] ==  'Pass') & (report['RunQC'] ==  'Pass')]
        
        # Run Number
        global run_number
        run_number = report['RunNumber'].unique()[0]
        
        # Adds sample name to respective list based of result and colours by SampleQC
        def process_sabin_items(ddns_classification, list_widget, classification_name):
            sabin_data = report.loc[report['DDNSclassification'].isin(ddns_classification)]

            if sabin_data['sample'].empty:
                list_widget.addItem(f'No {classification_name} present in Run')
                return
            for sabin in sabin_data['sample']:
                list_widget.addItem(sabin)
                items = list_widget.findItems(sabin, Qt.MatchExactly)

                if items:
                    item = items[0]
                    if sabin in ddns_pass['sample'].values:
                        #print(f'{sabin} has passed')
                        item.setBackground(QBrush(QColor('lightgreen')))
                    else:
                        #print(f'{sabin} has failed')
                        item.setBackground(QBrush(QColor('#f14141')))
                    
                    # initialise color tracking
                    item_color_changes[sabin] = item.background().color().name()

        
        # Global Dictionary to track item color changes
        global item_color_changes
        item_color_changes = {} 
                
        # # User colour change function
        # def change_item_color(item):
        #     current_color = item.background().color()
        #     base_text = item.text().split(' (')[0] # removes added text

        #     # Toggles between green and red on double-click and updates text
        #     # To Fail
        #     if current_color == QColor('lightgreen'):
        #         item.setBackground(QBrush(QColor('#f14141')))
        #         item.setText(f"{base_text} (To Fail)")  
        #         item_color_changes[item.text()] = '#f14141'
        #     # To Pass
        #     else:
        #         item.setBackground(QBrush(QColor('lightgreen')))
        #         item.setText(f"{base_text} (To Pass)")
        #         item_color_changes[item.text()] = 'lightgreen'

        #     print(f"Item '{base_text}' color changed to {item_color_changes[item.text()]}")
                            
        # Process each POLIO type
        process_sabin_items(['SABIN1', 'VDPV1'], self.s1_list, 'Sabin 1')
        process_sabin_items(['SABIN2', 'VDPV2'], self.s2_list, 'Sabin 2')
        process_sabin_items(['SABIN3', 'VDPV3'], self.s3_list, 'Sabin 3')

        # enable gen button after processing
        self.btn_save.setEnabled(True)

        # double-click signal to the color changing function
        # self.s1_list.itemDoubleClicked.connect(change_item_color)
        # self.s2_list.itemDoubleClicked.connect(change_item_color)
        # self.s3_list.itemDoubleClicked.connect(change_item_color)
        
        #print(item_color_changes)
        

    def save_fasta(self):

        # Iterating through FASTA files and extracting positive records
        output_fasta = fasta.split('/')[-1].replace('.fasta','_QCPass.fasta')
        db_fasta     = fasta.split('/')[-1].replace('.fasta','_DB.fasta')
        
        with open(f'{destination_path}/{output_fasta}', 'w') as output_qc, open(f'{destination_path}/{db_fasta}', 'w') as output_db:
            
            # Writer Class to ensure all sequences are wrapped at 60 nt and writes FASTAs to output file
            writer_qc = FastaWriter(output_qc, wrap=60)
            writer_qc.write_header()
            
            writer_db = FastaWriter(output_db, wrap=60)
            writer_db.write_header()

            for sequence in SeqIO.parse(fasta, 'fasta'):
                
                # extracts header sample name, and corrects if -DDNS is present
                split_header = sequence.id.rsplit('|')
                sample_id = split_header[0].replace('-DDNS','')
                
                #split_header[0] = sample_id
                
                # New Header elements
                # Extract barcode and ddns_group using regex
                barcode_match = re.search(r'barcode=(\S+)', sequence.description)
                ddns_group_match = re.search(r'ddns_group=(\S+)', sequence.description)

                if barcode_match and ddns_group_match:
                    barcode    = barcode_match.group(1)
                    ddns_group = ddns_group_match.group(1)
                else:
                    print(f"Error: barcode or ddns_group not found in sequence {sequence.id}")
                    continue  # Skip this sequence if elements are missing
                
                method = 'DDNS' if self.radio_ddns.isChecked() else 'Minion'
    
                # set new header
                new_header = [sample_id, run_number, barcode, ddns_group, 'QCPass', method]
                #print(new_header)
                fixed_header = '|'.join(new_header)
                
                # removes Nonpolio seqs
                if ddns_group == 'NonPolioEV':
                    continue
                
                # removes id that don't match sample names, removing differently named controls
                if not re.match(r'^[A-Z]{3}-\d{2}-\d{3,5}$', sample_id):
                    continue
                
                #print(ddns_group, barcode, sample_id)
                
                # Writes straight to database version with new header
                new_id_seq = SeqRecord(Seq(sequence.seq), id=fixed_header, description=sequence.description.split('|')[-1].strip())
                writer_db.write_record(new_id_seq)
                
                # checks if present in possible sample
                if sample_id in item_color_changes.keys():
                    # checks if green, aka Pass
                    if item_color_changes[sample_id] == '#90ee90':
            
                        # Save Final Record
                        writer_qc.write_record(new_id_seq)
                        #print(new_id_seq)
        

        QMessageBox.information(self, 'Success                                                                                                                               ', 
                                f'Samples that have passed QC have been saved to {output_fasta}.\n\nAll polio sequences for phylogenetic database have been saved to {db_fasta}.')
        pass
            
    def clear_list(self):
        self.listbox_view.clear()
        self.s1_list.clear()
        self.s2_list.clear()
        self.s3_list.clear()
        self.btn_save.setEnabled(False)

        # self.destination_entry.clear()
        # bg_label.setText('Drop CSVs here')
        # bg_label.show()
        
        
if __name__ == '__main__':
    app = QApplication(sys.argv)
    prog = App()
    prog.show()

    sys.exit(app.exec())
