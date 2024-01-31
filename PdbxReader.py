##
# File:  PdbxReader.py
# Original: 02-Feb-2009   jdw
# Date:  2012-01-09  Jdw  Adapted from PdbxParser
#
# Updates:
#
#   31-Jan-2024    Adrià Cereto Massagué merged PdbxContainers.py for convenience
#   12-Jan-2024    Adrià Cereto Massagué Made PEP-8 compliant, which fixes Jython compatibility
#   23-Mar-2011   jdw Added method to rename attributes in category containers.
#   05-Apr-2011   jdw Change cif writer to select double quoting as preferred
#                     quoting style where possible.
#   16-Jan-2012   jdw Create base class for DataCategory class
#   22-Mar-2012   jdw when append attributes to existing categories update
#                     existing rows with placeholder null values.
#    2-Sep-2012   jdw add option to avoid embedded quoting that might
#                     confuse simple parsers.
#   28-Jun-2013   jdw export remove method
#   29-Jun-2013   jdw export remove row method
# 2012-01-09 - (jdw) Separate reader and writer classes.
#
# 2012-09-02 - (jdw)  Revise tokenizer to better handle embedded quoting.
##
"""

A collection of container classes supporting the PDBx/mmCIF storage model
and PDBx/mmCIF dictionary and data file parser.

A base container class is defined which supports common features of
data and definition containers.   PDBx data files are organized in
sections called data blocks which are mapped to data containers.
PDBx dictionaries contain definition sections and data sections
which are mapped to definition and data containes respectively.

Data in both PDBx data files and dictionaries are organized in
data categories. In the PDBx syntax individual items or data
identified by labels of the form '_categoryName.attributeName'.
The terms category and attribute in PDBx jargon are analogous
table and column in relational data model, or class and attribute
in an object oriented data model.

The DataCategory class provides base storage container for instance
data and definition meta data.

Acknowledgements:

 The tokenizer used in this module is modeled after the clever parser design
 used in the PyMMLIB package.
 
 PyMMLib Development Group
 Authors: Ethan Merritt: merritt@u.washington.ed  & Jay Painter: jay.painter@gmail.com
 See:  http://pymmlib.sourceforge.net/

"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"

import re
import sys
import traceback


class CifName(object):
    ''' Class of utilities for CIF-style data names -
    '''

    def __init__(self):
        pass

    @staticmethod
    def categoryPart(name):
        tname = ""
        if name.startswith("_"):
            tname = name[1:]
        else:
            tname = name

        i = tname.find(".")
        if i == -1:
            return tname
        else:
            return tname[:i]

    @staticmethod
    def attributePart(name):
        i = name.find(".")
        if i == -1:
            return None
        else:
            return name[i + 1:]


class ContainerBase(object):
    ''' Container base class for data and definition objects.
    '''

    def __init__(self, name):
        # The enclosing scope of the data container (e.g. data_/save_)
        self.__name = name
        # List of category names within this container -
        self.__objNameList = []
        # dictionary of DataCategory objects keyed by category name.
        self.__objCatalog = {}
        self.__type = None

    def getType(self):
        return self.__type

    def setType(self, type):
        self.__type = type

    def getName(self):
        return self.__name

    def setName(self, name):
        self.__name = name

    def exists(self, name):
        if name in self.__objCatalog:
            return True
        else:
            return False

    def getObj(self, name):
        if name in self.__objCatalog:
            return self.__objCatalog[name]
        else:
            return None

    def getObjNameList(self):
        return self.__objNameList

    def append(self, obj):
        """ Add the input object to the current object catalog. An existing object
            of the same name will be overwritten.
        """
        if obj.getName() is not None:
            if obj.getName() not in self.__objCatalog:
                # self.__objNameList is keeping track of object order here --
                self.__objNameList.append(obj.getName())
            self.__objCatalog[obj.getName()] = obj

    def replace(self, obj):
        """ Replace an existing object with the input object
        """
        if ((obj.getName() is not None) and (
                obj.getName() in self.__objCatalog)):
            self.__objCatalog[obj.getName()] = obj

    def printIt(self, fh=sys.stdout, type="brief"):
        fh.write("+ %s container: %30s contains %4d categories\n" %
                 (self.getType(), self.getName(), len(self.__objNameList)))
        for nm in self.__objNameList:
            fh.write("--------------------------------------------\n")
            fh.write("Data category: %s\n" % nm)
            if type == 'brief':
                self.__objCatalog[nm].printIt(fh)
            else:
                self.__objCatalog[nm].dumpIt(fh)

    def rename(self, curName, newName):
        """ Change the name of an object in place -
        """
        try:
            i = self.__objNameList.index(curName)
            self.__objNameList[i] = newName
            self.__objCatalog[newName] = self.__objCatalog[curName]
            self.__objCatalog[newName].setName(newName)
            return True
        except BaseException:
            return False

    def remove(self, curName):
        """ Revmove object by name.  Return True on success or False otherwise.
        """
        try:
            if curName in self.__objCatalog:
                del self.__objCatalog[curName]
                i = self.__objNameList.index(curName)
                del self.__objNameList[i]
                return True
            else:
                return False
        except BaseException:
            pass

        return False


class DefinitionContainer(ContainerBase):
    def __init__(self, name):
        super(DefinitionContainer, self).__init__(name)
        self.setType('definition')

    def isCategory(self):
        if self.exists('category'):
            return True
        return False

    def isAttribute(self):
        if self.exists('item'):
            return True
        return False

    def printIt(self, fh=sys.stdout, type="brief"):
        fh.write("Definition container: %30s contains %4d categories\n" %
                 (self.getName(), len(self.getObjNameList())))
        if self.isCategory():
            fh.write("Definition type: category\n")
        elif self.isAttribute():
            fh.write("Definition type: item\n")
        else:
            fh.write("Definition type: undefined\n")

        for nm in self.getObjNameList():
            fh.write("--------------------------------------------\n")
            fh.write("Definition category: %s\n" % nm)
            if type == 'brief':
                self.getObj(nm).printIt(fh)
            else:
                self.getObj(nm).dumpId(fh)


class DataContainer(ContainerBase):
    ''' Container class for DataCategory objects.
    '''

    def __init__(self, name):
        super(DataContainer, self).__init__(name)
        self.setType('data')
        self.__globalFlag = False

    def invokeDataBlockMethod(self, type, method, db):
        self.__currentRow = 1
        exec(method.getInline())

    def setGlobal(self):
        self.__globalFlag = True

    def getGlobal(self):
        return self.__globalFlag


class DataCategoryBase(object):
    """ Base object definition for a data category -
    """

    def __init__(self, name, attributeNameList=None, rowList=None):
        self._name = name
        #
        if rowList is not None:
            self._rowList = rowList
        else:
            self._rowList = []

        if attributeNameList is not None:
            self._attributeNameList = attributeNameList
        else:
            self._attributeNameList = []
        #
        # Derived class data -
        #
        self._catalog = {}
        self._numAttributes = 0
        #
        self.__setup()

    def __setup(self):
        self._numAttributes = len(self._attributeNameList)
        self._catalog = {}
        for attributeName in self._attributeNameList:
            attributeNameLC = attributeName.lower()
            self._catalog[attributeNameLC] = attributeName
    #

    def setRowList(self, rowList):
        self._rowList = rowList

    def setAttributeNameList(self, attributeNameList):
        self._attributeNameList = attributeNameList
        self.__setup()

    def setName(self, name):
        self._name = name

    def get(self):
        return (self._name, self._attributeNameList, self._rowList)


class DataCategory(DataCategoryBase):
    """  Methods for creating, accessing, and formatting PDBx cif data categories.
    """

    def __init__(self, name, attributeNameList=None, rowList=None):
        super(DataCategory, self).__init__(name, attributeNameList, rowList)
        #
        self.__lfh = sys.stdout

        self.__currentRowIndex = 0
        self.__currentAttribute = None
        #
        self.__avoidEmbeddedQuoting = False
        #
        # --------------------------------------------------------------------
        # any whitespace
        self.__wsRe = re.compile(r"\s")
        self.__wsAndQuotesRe = re.compile(r"[\s'\"]")
        # any newline or carriage control
        self.__nlRe = re.compile(r"[\n\r]")
        #
        # single quote
        self.__sqRe = re.compile(r"[']")
        #
        self.__sqWsRe = re.compile(r"('\s)|(\s')")

        # double quote
        self.__dqRe = re.compile(r'["]')
        self.__dqWsRe = re.compile(r'("\s)|(\s")')
        #
        self.__intRe = re.compile(r'^[0-9]+$')
        self.__floatRe = re.compile(
            r'^-?(([0-9]+)[.]?|([0-9]*[.][0-9]+))([(][0-9]+[)])?([eE][+-]?[0-9]+)?$')
        #
        self.__dataTypeList = [
            'DT_NULL_VALUE',
            'DT_INTEGER',
            'DT_FLOAT',
            'DT_UNQUOTED_STRING',
            'DT_ITEM_NAME',
            'DT_DOUBLE_QUOTED_STRING',
            'DT_SINGLE_QUOTED_STRING',
            'DT_MULTI_LINE_STRING']
        self.__formatTypeList = [
            'FT_NULL_VALUE',
            'FT_NUMBER',
            'FT_NUMBER',
            'FT_UNQUOTED_STRING',
            'FT_QUOTED_STRING',
            'FT_QUOTED_STRING',
            'FT_QUOTED_STRING',
            'FT_MULTI_LINE_STRING']
        #

    def __getitem__(self, x):
        """  Implements list-type functionality -
             Implements op[x] for some special cases -
                x=integer - returns the row in category (normal list behavior)
                x=string  - returns the value of attribute 'x' in first row.
        """
        if isinstance(x, int):
            # return self._rowList.__getitem__(x)
            return self._rowList[x]

        elif isinstance(x, str):
            try:
                # return self._rowList[0][x]
                ii = self.getAttributeIndex(x)
                return self._rowList[0][ii]
            except (IndexError, KeyError):
                raise KeyError
        raise TypeError(x)

    def getCurrentAttribute(self):
        return self.__currentAttribute

    def getRowIndex(self):
        return self.__currentRowIndex

    def getRowList(self):
        return self._rowList

    def getRowCount(self):
        return (len(self._rowList))

    def getRow(self, index):
        try:
            return self._rowList[index]
        except BaseException:
            return []

    def removeRow(self, index):
        try:
            if ((index >= 0) and (index < len(self._rowList))):
                del self._rowList[index]
                if self.__currentRowIndex >= len(self._rowList):
                    self.__currentRowIndex = len(self._rowList) - 1
                return True
            else:
                pass
        except BaseException:
            pass

        return False

    def getFullRow(self, index):
        """ Return a full row based on the length of the the attribute list.
        """
        try:
            if (len(self._rowList[index]) < self._numAttributes):
                for ii in range(self._numAttributes -
                                len(self._rowList[index])):
                    self._rowList[index].append('?')
            return self._rowList[index]
        except BaseException:
            return ['?' for ii in range(self._numAttributes)]

    def getName(self):
        return self._name

    def getAttributeList(self):
        return self._attributeNameList

    def getAttributeCount(self):
        return len(self._attributeNameList)

    def getAttributeListWithOrder(self):
        oL = []
        for ii, att in enumerate(self._attributeNameList):
            oL.append((att, ii))
        return oL

    def getAttributeIndex(self, attributeName):
        try:
            return self._attributeNameList.index(attributeName)
        except BaseException:
            return -1

    def hasAttribute(self, attributeName):
        return attributeName in self._attributeNameList

    def getIndex(self, attributeName):
        try:
            return self._attributeNameList.index(attributeName)
        except BaseException:
            return -1

    def getItemNameList(self):
        itemNameList = []
        for att in self._attributeNameList:
            itemNameList.append("_" + self._name + "." + att)
        return itemNameList

    def append(self, row):
        # self.__lfh.write("PdbxContainer(append) category %s row %r\n"  % (self._name,row))
        self._rowList.append(row)

    def appendAttribute(self, attributeName):
        attributeNameLC = attributeName.lower()
        if attributeNameLC in self._catalog:
            i = self._attributeNameList.index(self._catalog[attributeNameLC])
            self._attributeNameList[i] = attributeName
            self._catalog[attributeNameLC] = attributeName
            # self.__lfh.write("Appending existing attribute %s\n" % attributeName)
        else:
            # self.__lfh.write("Appending existing attribute %s\n" % attributeName)
            self._attributeNameList.append(attributeName)
            self._catalog[attributeNameLC] = attributeName
            #
        self._numAttributes = len(self._attributeNameList)

    def appendAttributeExtendRows(self, attributeName):
        attributeNameLC = attributeName.lower()
        if attributeNameLC in self._catalog:
            i = self._attributeNameList.index(self._catalog[attributeNameLC])
            self._attributeNameList[i] = attributeName
            self._catalog[attributeNameLC] = attributeName
            self.__lfh.write(
                "Appending existing attribute %s\n" %
                attributeName)
        else:
            self._attributeNameList.append(attributeName)
            self._catalog[attributeNameLC] = attributeName
            # add a placeholder to any existing rows for the new attribute.
            if (len(self._rowList) > 0):
                for row in self._rowList:
                    row.append("?")
            #
        self._numAttributes = len(self._attributeNameList)

    def getValue(self, attributeName=None, rowIndex=None):
        if attributeName is None:
            attribute = self.__currentAttribute
        else:
            attribute = attributeName
        if rowIndex is None:
            rowI = self.__currentRowIndex
        else:
            rowI = rowIndex

        if isinstance(attribute, str) and isinstance(rowI, int):
            try:
                return self._rowList[rowI][self._attributeNameList.index(
                    attribute)]
            except (IndexError):
                raise IndexError
        raise IndexError(attribute)

    def setValue(self, value, attributeName=None, rowIndex=None):
        if attributeName is None:
            attribute = self.__currentAttribute
        else:
            attribute = attributeName

        if rowIndex is None:
            rowI = self.__currentRowIndex
        else:
            rowI = rowIndex

        if isinstance(attribute, str) and isinstance(rowI, int):
            try:
                # if row index is out of range - add the rows -
                for ii in range(rowI + 1 - len(self._rowList)):
                    self._rowList.append(self.__emptyRow())
                # self._rowList[rowI][attribute]=value
                ll = len(self._rowList[rowI])
                ind = self._attributeNameList.index(attribute)

                # extend the list if needed -
                if (ind >= ll):
                    self._rowList[rowI].extend(
                        [None for ii in xrange(2 * ind - ll)])
                self._rowList[rowI][ind] = value
            except (IndexError):
                self.__lfh.write(
                    "DataCategory(setvalue) index error category %s attribute %s index %d value %r\n" %
                    (self._name, attribute, rowI, value))
                traceback.print_exc(file=self.__lfh)
                # raise IndexError
            except (ValueError):
                self.__lfh.write(
                    "DataCategory(setvalue) value error category %s attribute %s index %d value %r\n" %
                    (self._name, attribute, rowI, value))
                traceback.print_exc(file=self.__lfh)
                # raise ValueError

    def __emptyRow(self):
        return [None for ii in range(len(self._attributeNameList))]

    def replaceValue(self, oldValue, newValue, attributeName):
        numReplace = 0
        if attributeName not in self._attributeNameList:
            return numReplace
        ind = self._attributeNameList.index(attributeName)
        for row in self._rowList:
            if row[ind] == oldValue:
                row[ind] = newValue
                numReplace += 1
        return numReplace

    def replaceSubstring(self, oldValue, newValue, attributeName):
        ok = False
        if attributeName not in self._attributeNameList:
            return ok
        ind = self._attributeNameList.index(attributeName)
        for row in self._rowList:
            val = row[ind]
            row[ind] = val.replace(oldValue, newValue)
            if val != row[ind]:
                ok = True
        return ok

    def invokeAttributeMethod(self, attributeName, type, method, db):
        self.__currentRowIndex = 0
        self.__currentAttribute = attributeName
        self.appendAttribute(attributeName)
        currentRowIndex = self.__currentRowIndex
        #
        ind = self._attributeNameList.index(attributeName)
        if len(self._rowList) == 0:
            row = [None for ii in xrange(len(self._attributeNameList) * 2)]
            row[ind] = None
            self._rowList.append(row)

        for row in self._rowList:
            ll = len(row)
            if (ind >= ll):
                row.extend([None for ii in xrange(2 * ind - ll)])
                row[ind] = None
            exec(method.getInline())
            self.__currentRowIndex += 1
            currentRowIndex = self.__currentRowIndex

    def invokeCategoryMethod(self, type, method, db):
        self.__currentRowIndex = 0
        exec(method.getInline())

    def getAttributeLengthMaximumList(self):
        mList = [0 for i in len(self._attributeNameList)]
        for row in self._rowList:
            for indx, val in enumerate(row):
                mList[indx] = max(mList[indx], len(val))
        return mList

    def renameAttribute(self, curAttributeName, newAttributeName):
        """ Change the name of an attribute in place -
        """
        try:
            i = self._attributeNameList.index(curAttributeName)
            self._attributeNameList[i] = newAttributeName
            del self._catalog[curAttributeName.lower()]
            self._catalog[newAttributeName.lower()] = newAttributeName
            return True
        except BaseException:
            return False

    def printIt(self, fh=sys.stdout):
        fh.write("--------------------------------------------\n")
        fh.write("  Category: %s attribute list length: %d\n" %
                 (self._name, len(self._attributeNameList)))
        for at in self._attributeNameList:
            fh.write("  Category: %s attribute: %s\n" % (self._name, at))

        fh.write("  Row value list length: %d\n" % len(self._rowList))
        #
        for row in self._rowList[:2]:
            #
            if len(row) == len(self._attributeNameList):
                for ii, v in enumerate(row):
                    fh.write("        %30s: %s ...\n" %
                             (self._attributeNameList[ii], str(v)[:30]))
            else:
                fh.write(
                    "+WARNING - %s data length %d attribute name length %s mismatched\n" %
                    (self._name, len(row), len(
                        self._attributeNameList)))

    def dumpIt(self, fh=sys.stdout):
        fh.write("--------------------------------------------\n")
        fh.write("  Category: %s attribute list length: %d\n" %
                 (self._name, len(self._attributeNameList)))
        for at in self._attributeNameList:
            fh.write("  Category: %s attribute: %s\n" % (self._name, at))

        fh.write("  Value list length: %d\n" % len(self._rowList))
        for row in self._rowList:
            for ii, v in enumerate(row):
                fh.write(
                    "        %30s: %s\n" %
                    (self._attributeNameList[ii], v))

    def __formatPdbx(self, inp):
        """ Format input data following PDBx quoting rules -
        """
        try:
            if (inp is None):
                return ("?", 'DT_NULL_VALUE')

            # pure numerical values are returned as unquoted strings
            if (isinstance(inp, int) or self.__intRe.search(str(inp))):
                return ([str(inp)], 'DT_INTEGER')

            if (isinstance(inp, float) or self.__floatRe.search(str(inp))):
                return ([str(inp)], 'DT_FLOAT')

            # null value handling -

            if (inp == "." or inp == "?"):
                return ([inp], 'DT_NULL_VALUE')

            if (inp == ""):
                return (["."], 'DT_NULL_VALUE')

            # Contains white space or quotes ?
            if not self.__wsAndQuotesRe.search(inp):
                if inp.startswith("_"):
                    return (self.__doubleQuotedList(inp), 'DT_ITEM_NAME')
                else:
                    return ([str(inp)], 'DT_UNQUOTED_STRING')
            else:
                if self.__nlRe.search(inp):
                    return (
                        self.__semiColonQuotedList(inp),
                        'DT_MULTI_LINE_STRING')
                else:
                    if (self.__avoidEmbeddedQuoting):
                        # change priority to choose double quoting where
                        # possible.
                        if not self.__dqRe.search(
                                inp) and not self.__sqWsRe.search(inp):
                            return (
                                self.__doubleQuotedList(inp),
                                'DT_DOUBLE_QUOTED_STRING')
                        elif not self.__sqRe.search(inp) and not self.__dqWsRe.search(inp):
                            return (
                                self.__singleQuotedList(inp),
                                'DT_SINGLE_QUOTED_STRING')
                        else:
                            return (
                                self.__semiColonQuotedList(inp),
                                'DT_MULTI_LINE_STRING')
                    else:
                        # change priority to choose double quoting where
                        # possible.
                        if not self.__dqRe.search(inp):
                            return (
                                self.__doubleQuotedList(inp),
                                'DT_DOUBLE_QUOTED_STRING')
                        elif not self.__sqRe.search(inp):
                            return (
                                self.__singleQuotedList(inp),
                                'DT_SINGLE_QUOTED_STRING')
                        else:
                            return (
                                self.__semiColonQuotedList(inp),
                                'DT_MULTI_LINE_STRING')

        except BaseException:
            traceback.print_exc(file=self.__lfh)

    def __dataTypePdbx(self, inp):
        """ Detect the PDBx data type -
        """
        if (inp is None):
            return ('DT_NULL_VALUE')

        # pure numerical values are returned as unquoted strings
        if isinstance(inp, int) or self.__intRe.search(str(inp)):
            return ('DT_INTEGER')

        if isinstance(inp, float) or self.__floatRe.search(str(inp)):
            return ('DT_FLOAT')

        # null value handling -

        if (inp == "." or inp == "?"):
            return ('DT_NULL_VALUE')

        if (inp == ""):
            return ('DT_NULL_VALUE')

        # Contains white space or quotes ?
        if not self.__wsAndQuotesRe.search(inp):
            if inp.startswith("_"):
                return ('DT_ITEM_NAME')
            else:
                return ('DT_UNQUOTED_STRING')
        else:
            if self.__nlRe.search(inp):
                return ('DT_MULTI_LINE_STRING')
            else:
                if (self.__avoidEmbeddedQuoting):
                    if not self.__sqRe.search(
                            inp) and not self.__dqWsRe.search(inp):
                        return ('DT_DOUBLE_QUOTED_STRING')
                    elif not self.__dqRe.search(inp) and not self.__sqWsRe.search(inp):
                        return ('DT_SINGLE_QUOTED_STRING')
                    else:
                        return ('DT_MULTI_LINE_STRING')
                else:
                    if not self.__sqRe.search(inp):
                        return ('DT_DOUBLE_QUOTED_STRING')
                    elif not self.__dqRe.search(inp):
                        return ('DT_SINGLE_QUOTED_STRING')
                    else:
                        return ('DT_MULTI_LINE_STRING')

    def __singleQuotedList(self, inp):
        l = []
        l.append("'")
        l.append(inp)
        l.append("'")
        return (l)

    def __doubleQuotedList(self, inp):
        l = []
        l.append('"')
        l.append(inp)
        l.append('"')
        return (l)

    def __semiColonQuotedList(self, inp):
        l = []
        l.append("\n")
        if inp[-1] == '\n':
            l.append(";")
            l.append(inp)
            l.append(";")
            l.append("\n")
        else:
            l.append(";")
            l.append(inp)
            l.append("\n")
            l.append(";")
            l.append("\n")

        return (l)

    def getValueFormatted(self, attributeName=None, rowIndex=None):
        if attributeName is None:
            attribute = self.__currentAttribute
        else:
            attribute = attributeName

        if rowIndex is None:
            rowI = self.__currentRowIndex
        else:
            rowI = rowIndex

        if isinstance(attribute, str) and isinstance(rowI, int):
            try:
                list, type = self.__formatPdbx(
                    self._rowList[rowI][self._attributeNameList.index(attribute)])
                return "".join(list)
            except (IndexError):
                self.__lfh.write(
                    "attributeName %s rowI %r rowdata %r\n" %
                    (attributeName, rowI, self._rowList[rowI]))
                raise IndexError
        raise TypeError(attribute)

    def getValueFormattedByIndex(self, attributeIndex, rowIndex):
        try:
            list, type = self.__formatPdbx(
                self._rowList[rowIndex][attributeIndex])
            return "".join(list)
        except (IndexError):
            raise IndexError

    def getAttributeValueMaxLengthList(self, steps=1):
        mList = [0 for i in range(len(self._attributeNameList))]
        for row in self._rowList[::steps]:
            for indx in range(len(self._attributeNameList)):
                val = row[indx]
                mList[indx] = max(mList[indx], len(str(val)))
        return mList

    def getFormatTypeList(self, steps=1):
        try:
            curDataTypeList = ['DT_NULL_VALUE' for i in range(
                len(self._attributeNameList))]
            for row in self._rowList[::steps]:
                for indx in range(len(self._attributeNameList)):
                    val = row[indx]
                    # print "index ",indx," val ",val
                    dType = self.__dataTypePdbx(val)
                    dIndx = self.__dataTypeList.index(dType)
                    # print "d type", dType, " d type index ",dIndx

                    cType = curDataTypeList[indx]
                    cIndx = self.__dataTypeList.index(cType)
                    cIndx = max(cIndx, dIndx)
                    curDataTypeList[indx] = self.__dataTypeList[cIndx]

            # Map the format types to the data types
            curFormatTypeList = []
            for dt in curDataTypeList:
                ii = self.__dataTypeList.index(dt)
                curFormatTypeList.append(self.__formatTypeList[ii])
        except BaseException:
            self.__lfh.write(
                "PdbxDataCategory(getFormatTypeList) ++Index error at index %d in row %r\n" %
                (indx, row))

        return curFormatTypeList, curDataTypeList

    def getFormatTypeListX(self):
        curDataTypeList = ['DT_NULL_VALUE' for i in range(
            len(self._attributeNameList))]
        for row in self._rowList:
            for indx in range(len(self._attributeNameList)):
                val = row[indx]
                # print "index ",indx," val ",val
                dType = self.__dataTypePdbx(val)
                dIndx = self.__dataTypeList.index(dType)
                # print "d type", dType, " d type index ",dIndx

                cType = curDataTypeList[indx]
                cIndx = self.__dataTypeList.index(cType)
                cIndx = max(cIndx, dIndx)
                curDataTypeList[indx] = self.__dataTypeList[cIndx]

        # Map the format types to the data types
        curFormatTypeList = []
        for dt in curDataTypeList:
            ii = self.__dataTypeList.index(dt)
            curFormatTypeList.append(self.__formatTypeList[ii])
        return curFormatTypeList, curDataTypeList


class PdbxError(Exception):
    """ Class for catch general errors 
    """
    pass

class SyntaxError(Exception):
    """ Class for catching syntax errors 
    """
    def __init__(self, lineNumber, text):
        Exception.__init__(self)
        self.lineNumber = lineNumber
        self.text = text

    def __str__(self):
        return "%%ERROR - [at line: %d] %s" % (self.lineNumber, self.text)



class PdbxReader(object):
    """ PDBx reader for data files and dictionaries.
    
    """
    def __init__(self,ifh):
        """  ifh - input file handle returned by open()
        """
        # 
        self.__curLineNumber = 0        
        self.__ifh=ifh
        self.__stateDict={"data":   "ST_DATA_CONTAINER",
                          "loop":   "ST_TABLE",
                          "global": "ST_GLOBAL_CONTAINER",
                          "save":   "ST_DEFINITION",
                          "stop":   "ST_STOP"}
        
    def read(self, containerList):
        """
        Appends to the input list of definition and data containers.
        
        """
        self.__curLineNumber = 0
        try:
            self.__parser(self.__tokenizer(self.__ifh), containerList)
        except StopIteration:
            pass
        else:
            raise PdbxError()

    def __syntaxError(self, errText):
        raise SyntaxError(self.__curLineNumber, errText)

    def __getContainerName(self,inWord):
        """ Returns the name of the data_ or save_ container
        """
        return str(inWord[5:]).strip()
    
    def __getState(self, inWord):
        """Identifies reserved syntax elements and assigns an associated state.  

           Returns: (reserved word, state)
           where - 
              reserved word -  is one of CIF syntax elements:
                               data_, loop_, global_, save_, stop_
              state - the parser state required to process this next section.
        """
        i = inWord.find("_")
        if i == -1:
            return None,"ST_UNKNOWN"

        try:
            rWord=inWord[:i].lower()            
            return rWord, self.__stateDict[rWord]
        except:
            return None,"ST_UNKNOWN"
        
    def __parser(self, tokenizer, containerList):
        """ Parser for PDBx data files and dictionaries.

            Input - tokenizer() reentrant method recognizing data item names (_category.attribute)
                    quoted strings (single, double and multi-line semi-colon delimited), and unquoted
                    strings.

                    containerList -  list-type container for data and definition objects parsed from
                                     from the input file.

            Return:
                    containerList - is appended with data and definition objects - 
        """
        # Working container - data or definition
        curContainer = None
        #
        # Working category container 
        categoryIndex = {}
        curCategory = None
        #
        curRow = None
        state =  None

        # Find the first reserved word and begin capturing data.
        #
        while True:
            curCatName, curAttName, curQuotedString, curWord = next(tokenizer)
            if curWord is None:
                continue
            reservedWord, state  = self.__getState(curWord)
            if reservedWord is not None:
                break
        
        while True:
            #
            #  Set the current state  -
            #
            #  At this point in the processing cycle we are expecting a token containing
            #  either a '_category.attribute'  or a reserved word.  
            #
            if curCatName is not None:
                state = "ST_KEY_VALUE_PAIR"
            elif curWord is not None:
                reservedWord, state = self.__getState(curWord)
            else:
                self.__syntaxError("Miscellaneous syntax error")
                return            

            #
            # Process  _category.attribute  value assignments 
            #
            if state == "ST_KEY_VALUE_PAIR":
                try:
                    curCategory = categoryIndex[curCatName]
                except KeyError:
                    # A new category is encountered - create a container and add a row 
                    curCategory = categoryIndex[curCatName] = DataCategory(curCatName)

                    try:
                        curContainer.append(curCategory)
                    except AttributeError:
                        self.__syntaxError("Category cannot be added to  data_ block")
                        return

                    curRow = []                    
                    curCategory.append(curRow)
                else:
                    # Recover the existing row from the category
                    try:
                        curRow = curCategory[0] 
                    except IndexError:
                        self.__syntaxError("Internal index error accessing category data")
                        return

                # Check for duplicate attributes and add attribute to table.
                if curAttName in curCategory.getAttributeList():
                    self.__syntaxError("Duplicate attribute encountered in category")
                    return
                else:
                    curCategory.appendAttribute(curAttName)


                # Get the data for this attribute from the next token
                tCat, tAtt, curQuotedString, curWord = next(tokenizer)

                if tCat is not None or (curQuotedString is None and curWord is None):
                    self.__syntaxError("Missing data for item _%s.%s" % (curCatName,curAttName))

                if curWord is not None:
                    # 
                    # Validation check token for misplaced reserved words  -  
                    #
                    reservedWord, state  = self.__getState(curWord)
                    if reservedWord is not None:
                        self.__syntaxError("Unexpected reserved word: %s" % (reservedWord))

                    curRow.append(curWord)

                elif curQuotedString is not None:
                    curRow.append(curQuotedString)

                else:
                    self.__syntaxError("Missing value in item-value pair")

                curCatName, curAttName, curQuotedString, curWord = next(tokenizer)
                continue

            #
            # Process a loop_ declaration and associated data -
            #
            elif state == "ST_TABLE":

                # The category name in the next curCatName,curAttName pair
                #    defines the name of the category container.
                curCatName,curAttName,curQuotedString,curWord = next(tokenizer)

                if curCatName is None or curAttName is None:
                    self.__syntaxError("Unexpected token in loop_ declaration")
                    return

                # Check for a previous category declaration.
                if curCatName in categoryIndex:
                    self.__syntaxError("Duplicate category declaration in loop_")
                    return

                curCategory = DataCategory(curCatName)

                try:
                    curContainer.append(curCategory)
                except AttributeError:
                    self.__syntaxError("loop_ declaration outside of data_ block or save_ frame")
                    return

                curCategory.appendAttribute(curAttName)

                # Read the rest of the loop_ declaration 
                while True:
                    curCatName, curAttName, curQuotedString, curWord = next(tokenizer)
                    
                    if curCatName is None:
                        break

                    if curCatName != curCategory.getName():
                        self.__syntaxError("Changed category name in loop_ declaration")
                        return

                    curCategory.appendAttribute(curAttName)


                # If the next token is a 'word', check it for any reserved words - 
                if curWord is not None:
                    reservedWord, state  = self.__getState(curWord)
                    if reservedWord is not None:
                        if reservedWord == "stop":
                            return
                        else:
                            self.__syntaxError("Unexpected reserved word after loop declaration: %s" % (reservedWord))
                    
                # Read the table of data for this loop_ - 
                while True:
                    curRow = []                    
                    curCategory.append(curRow)

                    for tAtt in curCategory.getAttributeList():
                        if curWord is not None:
                            curRow.append(curWord)
                        elif curQuotedString is not None:
                            curRow.append(curQuotedString)

                        curCatName,curAttName,curQuotedString,curWord = next(tokenizer)

                    # loop_ data processing ends if - 

                    # A new _category.attribute is encountered
                    if curCatName is not None:
                        break

                    # A reserved word is encountered
                    if curWord is not None:
                        reservedWord, state = self.__getState(curWord)
                        if reservedWord is not None:
                            break
                        
                continue


            elif state == "ST_DEFINITION":
                # Ignore trailing unnamed saveframe delimiters e.g. 'save_'
                sName=self.__getContainerName(curWord)
                if (len(sName) > 0):
                    curContainer = DefinitionContainer(sName)
                    containerList.append(curContainer)
                    categoryIndex = {}
                    curCategory = None

                curCatName,curAttName,curQuotedString,curWord = next(tokenizer)

            elif state == "ST_DATA_CONTAINER":
                #
                dName=self.__getContainerName(curWord)
                if len(dName) == 0:
                    dName="unidentified"
                curContainer = DataContainer(dName)
                containerList.append(curContainer)
                categoryIndex = {}
                curCategory = None
                curCatName,curAttName,curQuotedString,curWord = next(tokenizer)

            elif state == "ST_STOP":
                return
            elif state == "ST_GLOBAL":
                curContainer = DataContainer("blank-global")
                curContainer.setGlobal()
                containerList.append(curContainer)
                categoryIndex = {}
                curCategory = None
                curCatName,curAttName,curQuotedString,curWord = next(tokenizer)

            elif state == "ST_UNKNOWN":
                self.__syntaxError("Unrecogized syntax element: " + str(curWord))
                return
                

    def __tokenizer(self, ifh):
        """ Tokenizer method for the mmCIF syntax file - 

            Each return/yield from this method returns information about
            the next token in the form of a tuple with the following structure.

            (category name, attribute name, quoted strings, words w/o quotes or white space)

            Differentiated the reqular expression to the better handle embedded quotes.

        """
        #
        # Regex definition for mmCIF syntax - semi-colon delimited strings are handled
        #                                     outside of this regex.
        mmcifRe = re.compile(
            r"(?:"

             "(?:_(.+?)[.](\S+))"               "|"  # _category.attribute

             "(?:['](.*?)(?:[']\s|[']$))"       "|"  # single quoted strings
             "(?:[\"](.*?)(?:[\"]\s|[\"]$))"    "|"  # double quoted strings             

             "(?:\s*#.*$)"                      "|"  # comments (dumped)

             "(\S+)"                                 # unquoted words

             ")")

        fileIter = iter(ifh)

        ## Tokenizer loop begins here ---
        while True:
            try:
                line = next(fileIter)
            except StopIteration:
                break
            self.__curLineNumber += 1

            # Dump comments
            if line.startswith("#"):
                continue
            
            # Gobble up the entire semi-colon/multi-line delimited string and
            #    and stuff this into the string slot in the return tuple
            #
            if line.startswith(";"):
                mlString = [line[1:]]
                while True:
                    line = next(fileIter)
                    self.__curLineNumber += 1
                    if line.startswith(";"):
                        break
                    mlString.append(line)

                # remove trailing new-line that is part of the \n; delimiter
                mlString[-1] = mlString[-1].rstrip()
                #
                yield (None, None, "".join(mlString), None)
                #
                # Need to process the remainder of the current line -
                line = line[1:]
                #continue

            # Apply regex to the current line consolidate the single/double
            # quoted within the quoted string category
            for it in mmcifRe.finditer(line):
                tgroups = it.groups()
                if tgroups != (None, None, None, None, None):
                    if tgroups[2] is not None:
                        qs = tgroups[2]
                    elif tgroups[3] is not None:
                        qs = tgroups[3]
                    else:
                        qs = None
                    groups = (tgroups[0],tgroups[1],qs,tgroups[4])
                    yield groups

    def __tokenizerOrg(self, ifh):
        """ Tokenizer method for the mmCIF syntax file - 

            Each return/yield from this method returns information about
            the next token in the form of a tuple with the following structure.

            (category name, attribute name, quoted strings, words w/o quotes or white space)

        """
        #
        # Regex definition for mmCIF syntax - semi-colon delimited strings are handled
        #                                     outside of this regex.
        mmcifRe = re.compile(
            r"(?:"

             "(?:_(.+?)[.](\S+))"               "|"  # _category.attribute

             "(?:['\"](.*?)(?:['\"]\s|['\"]$))" "|"  # quoted strings

             "(?:\s*#.*$)"                      "|"  # comments (dumped)

             "(\S+)"                                 # unquoted words

             ")")

        fileIter = iter(ifh)

        ## Tokenizer loop begins here ---
        while True:
            line = next(fileIter)
            self.__curLineNumber += 1

            # Dump comments
            if line.startswith("#"):
                continue
            
            # Gobble up the entire semi-colon/multi-line delimited string and
            #    and stuff this into the string slot in the return tuple
            #
            if line.startswith(";"):
                mlString = [line[1:]]
                while True:
                    line = next(fileIter)
                    self.__curLineNumber += 1
                    if line.startswith(";"):
                        break
                    mlString.append(line)

                # remove trailing new-line that is part of the \n; delimiter
                mlString[-1] = mlString[-1].rstrip()
                #
                yield (None, None, "".join(mlString), None)
                #
                # Need to process the remainder of the current line -
                line = line[1:]
                #continue

            ## Apply regex to the current line 
            for it in mmcifRe.finditer(line):
                groups = it.groups()
                if groups != (None, None, None, None):
                    yield groups
