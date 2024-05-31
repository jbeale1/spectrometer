# edit a set of parameters in one table
# https://learndataanalysis.org/create-a-pandas-dataframe-editor-with-pyqt5/

import sys
import pandas as pd
from PyQt5.QtWidgets import QApplication, QWidget, QTableWidget, QTableWidgetItem, QLineEdit, \
                            QPushButton, QItemDelegate, QVBoxLayout
from PyQt5.QtGui import QDoubleValidator

class FloatDelegate(QItemDelegate):
    def __init__(self, parent=None):
        super().__init__()

    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        editor.setValidator(QDoubleValidator())
        return editor

class TableWidget(QTableWidget):
    def __init__(self, df):
        super().__init__()
        self.df = df
        self.resize(380, 150)
        self.setStyleSheet('font-size: 14px;')

        # set table dimension
        nRows, nColumns = self.df.shape
        self.setColumnCount(nColumns)
        self.setRowCount(nRows)
        self.setHorizontalHeaderLabels(('Parameter', 'Value'))
        self.setColumnWidth(0, 150)
        self.setColumnWidth(1, 70)
        self.setItemDelegateForColumn(1, FloatDelegate())

        # data insertion
        for i in range(self.rowCount()):
            for j in range(self.columnCount()):
                self.setItem(i, j, QTableWidgetItem(str(self.df.iloc[i, j])))

        self.cellChanged[int, int].connect(self.updateDF)   

    def updateDF(self, row, column):
        text = self.item(row, column).text()
        self.df.iloc[row, column] = text

class ParamEditor(QWidget):

    def __init__(self, df):
        super().__init__()
        self.df = df

        self.resize(300, 250)

        mainLayout = QVBoxLayout()

        self.table = TableWidget(self.df)
        mainLayout.addWidget(self.table)

        button_print = QPushButton('Display DF')
        button_print.setStyleSheet('font-size: 14px')
        button_print.clicked.connect(self.print_DF_Values)
        mainLayout.addWidget(button_print)

        button_export = QPushButton('Save Param File')
        button_export.setStyleSheet('font-size: 14px')
        button_export.clicked.connect(self.export_to_csv)
        mainLayout.addWidget(button_export)     

        self.setLayout(mainLayout)
        
    def print_DF_Values(self):
        df = self.table.df
        print(df)
        print("scanMin = ", df.loc[['scanMin']].values[0][1] )
        print("scanMax = ", df.loc[['scanMax']].values[0][1] )

    def export_to_csv(self):
        self.table.df.to_csv('SpecParams.csv', index=False)
        print('Parameter file exported.')

if __name__ == '__main__':
    app = QApplication(sys.argv)


    paramList = ['scanMin','scanMax','pkStart','pkStop']
    data = {
        'Parameter': paramList,
        'Value': [360, 800, 385, 750]
    }
    df = pd.DataFrame(data)
    df.index = paramList # use names for row index, not just the default numbers

    demo = ParamEditor(df)
    demo.show()
    retval = app.exec_()

    print("Window has been closed, now dataframe is:
          ")
    print(df)
    sys.exit(retval)
    
