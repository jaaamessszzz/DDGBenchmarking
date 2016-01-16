#!/usr/bin/python2.4
# encoding: utf-8
"""
db_schema.py
SQLAlchemy representation

Created by Shane O'Connor 2015.
Copyright (c) 2015 Shane O'Connor. All rights reserved.
"""

import sys
import os
import inspect
import pprint

from sqlalchemy import Table, Column, Integer, ForeignKey
from sqlalchemy import inspect as sqlalchemy_inspect
from sqlalchemy.orm import relationship, backref
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.dialects.mysql import DOUBLE, TINYINT, LONGBLOB
from sqlalchemy.types import DateTime, Enum, Integer, TIMESTAMP, Text, Unicode, String


if __name__ == '__main__':
    sys.path.insert(0, '../../klab')

from klab.db.sqlalchemy_interface import MySQLSchemaConverter
from klab.fs.fsio import read_file
from klab import colortext


DeclarativeBase = declarative_base()


#############################
#                           #
#  Foundational constructs  #
#                           #
#############################



class AminoAcid(DeclarativeBase):
    __tablename__ = 'AminoAcid'

    Code = Column(String(1), nullable=False, primary_key=True)
    LongCode = Column(String(3), nullable=False)
    Name = Column(String(32), nullable=False)
    Polarity = Column(Enum('polar','non-polar','charged'), nullable=True)
    Aromaticity = Column(Enum('aromatic','aliphatic','neither'), nullable=True)
    Hydrophobicity_pH7 = Column(Enum('hydrophobic','hydrophilic'), nullable=True)
    SideChainAcidity = Column(Enum('acidic','basic','neutral'), nullable=True)
    pKa = Column(DOUBLE, nullable=True)
    AverageMass = Column(DOUBLE, nullable=True)
    Volume = Column(DOUBLE, nullable=True)
    Size = Column(Enum('small','large'), nullable=True)
    Tiny = Column(TINYINT(1), nullable=True, default=0)


class FileContent(DeclarativeBase):
    __tablename__ = 'FileContent'

    ID = Column(Integer, nullable=False, primary_key=True)
    Content = Column(LONGBLOB, nullable=False)
    MIMEType = Column(String(64), nullable=False)
    Filesize = Column(Integer, nullable=False)
    MD5HexDigest = Column(String(32), nullable=False)


######################################
#                                    #
#  Ligands and associated records    #
#                                    #
######################################


class Ligand(DeclarativeBase):
    __tablename__ = 'Ligand'

    ID = Column(Integer, nullable=False, primary_key=True)
    PDBCode = Column(String(3), nullable=False)
    LigandCode = Column(String(512), nullable=False)
    Formula = Column(String(256), nullable=False)
    MolecularWeight = Column(DOUBLE, nullable=False)
    LigandType = Column(String(256), nullable=False)
    Solubility = Column(String(256), nullable=True)
    CellPermeability = Column(Enum('Yes','No','Yes (probably)'), nullable=True)
    AssaysToDetermineConcentrationInCells = Column(String(256), nullable=True)
    ProductionInCells = Column(Enum('Yes','No'), nullable=True)
    ProductionInCellsNotes = Column(String(64), nullable=True)
    Diagram = Column(LONGBLOB, nullable=True)
    SimilarCompoundsDiagram = Column(LONGBLOB, nullable=True)
    InChI = Column(Text, nullable=False)
    InChIKey = Column(String(27), nullable=False)


class LigandDescriptor(DeclarativeBase):
    __tablename__ = 'LigandDescriptor'

    ID = Column(Integer, nullable=False, primary_key=True)
    LigandID = Column(Integer, ForeignKey('Ligand.ID'), nullable=False)
    Descriptor = Column(Text, nullable=False)
    DescriptorType = Column(String(128), nullable=False)
    Program = Column(String(64), nullable=False)
    Version = Column(String(16), nullable=False)


class LigandIdentifier(DeclarativeBase):
    __tablename__ = 'LigandIdentifier'

    ID = Column(Integer, nullable=False, primary_key=True)
    LigandID = Column(Integer, ForeignKey('Ligand.ID'), nullable=False)
    Identifier = Column(String(1024), nullable=False)
    IDType = Column(String(128), nullable=False)
    Program = Column(String(64), nullable=False)
    Version = Column(String(16), nullable=False)


class LigandPrice(DeclarativeBase):
    __tablename__ = 'LigandPrice'

    LigandID = Column(Integer, ForeignKey('Ligand.ID'), nullable=False, primary_key=True)
    PriceDate = Column(DateTime, nullable=False, primary_key=True)
    USDPricePerGram = Column(DOUBLE, nullable=False)
    PriceNote = Column(String(256), nullable=True)


class LigandReference(DeclarativeBase):
    __tablename__ = 'LigandReference'

    ID = Column(Integer, nullable=False, primary_key=True)
    LigandID = Column(Integer, ForeignKey('Ligand.ID'), nullable=False)
    PublicationID = Column(String(64), ForeignKey('Publication.ID'), nullable=True)
    Type = Column(Enum('Reference','Assay'), nullable=False, default=u'Reference')
    Notes = Column(Text, nullable=True)


class LigandSynonym(DeclarativeBase):
    __tablename__ = 'LigandSynonym'

    LigandID = Column(Integer, ForeignKey('Ligand.ID'), nullable=False, primary_key=True)
    Synonym = Column(String(256), nullable=False, primary_key=True)


######################################
#                                    #
#  Ions and associated records       #
#                                    #
######################################


class Ion(DeclarativeBase):
    __tablename__ = 'Ion'

    ID = Column(Integer, nullable=False, primary_key=True)
    PDBCode = Column(String(3), nullable=False)
    Formula = Column(String(256), nullable=False)
    Description = Column(String(256), nullable=True)


######################################
#                                    #
#  PDB files and associated records  #
#                                    #
######################################


class PDBFile(DeclarativeBase):
    __tablename__ = 'PDBFile'

    ID = Column(String(10), nullable=False, primary_key=True, default=u'')
    FileSource = Column(String(64), nullable=False, default=u'RCSB')
    Content = Column(Text, nullable=False)
    FASTA = Column(Text, nullable=False)
    Resolution = Column(DOUBLE, nullable=True)
    Techniques = Column(String(256), nullable=False)
    BFactorMean = Column(DOUBLE, nullable=True)
    BFactorDeviation = Column(DOUBLE, nullable=True)
    Publication = Column(String(64), ForeignKey('Publication.ID'), nullable=True)
    Transmembrane = Column(TINYINT(1), nullable=True)
    UserID = Column(String(64), nullable=True)
    Notes = Column(Text, nullable=True)
    DerivedFrom = Column(String(4), nullable=True)

    # Relationships
    residues = relationship('Publication', primaryjoin="PDBFile.Publication==Publication.ID")

    def __repr__(self):
        notes = ''
        resolution = 'unknown resolution'
        if self.Resolution:
            resolution = str(self.Resolution) + 'A'
        if self.Notes and self.Notes.strip():
            notes = self.Notes.strip()
            if notes[-1] != '.':
                notes += '.'
        return 'PDBFile: {0}, {1} ({2}). Source: {3}. B-Factors: {4} (mean), {5} (stddev). Transmembrane protein: {6}. {7}'.format(self.ID, resolution, self.Techniques, self.FileSource, self.BFactorMean, self.BFactorDeviation, self.Transmembrane, notes)


class PDBChain(DeclarativeBase):
    __tablename__ = 'PDBChain'

    PDBFileID = Column(String(10), ForeignKey('PDBFile.ID'), nullable=False, primary_key=True)
    Chain = Column(String(1), nullable=False, primary_key=True)
    MoleculeType = Column(Enum('Protein', 'DNA', 'RNA', 'Ligand', 'Protein skeleton', 'Heterogen', 'Solution', 'Unknown'), nullable=True)
    WildtypeProteinID = Column(String(18), nullable=True)
    FullProteinID = Column(String(18), nullable=True)
    SegmentProteinID = Column(String(18), nullable=True)
    WildtypeAlignedProteinID = Column(String(18), nullable=True)
    AcquiredProteinID = Column(String(18), nullable=True)
    Coordinates = Column(LONGBLOB, nullable=True)

    # Parent relationships
    pdb_file = relationship('PDBFile', primaryjoin="PDBChain.PDBFileID==PDBFile.ID")

    # Children relationships
    residues = relationship('PDBResidue', primaryjoin="and_(PDBResidue.PDBFileID==PDBChain.PDBFileID, PDBResidue.Chain==PDBChain.Chain)")

    def __init__(self, **kwargs):
        super(PDBChain, self).__init__(**kwargs)
        # do custom initialization here

    def __repr__(self):
        return 'PDBChain: {0}, {1} ({2})'.format(self.PDBFileID, self.Chain, self.MoleculeType)


class PDBMolecule(DeclarativeBase):
    __tablename__ = 'PDBMolecule'

    PDBFileID = Column(String(10), ForeignKey('PDBFile.ID'), nullable=False, primary_key=True)
    MoleculeID = Column(Integer, nullable=False, primary_key=True)
    Name = Column(String(256), nullable=False)
    Organism = Column(String(256), nullable=True)
    Fragment = Column(String(256), nullable=True)
    Synonym = Column(String(256), nullable=True)
    Engineered = Column(TINYINT(1), nullable=True)
    EC = Column(String(32), nullable=True)
    Mutation = Column(TINYINT(1), nullable=True)
    OtherDetails = Column(String(256), nullable=True)

    # Parent relationships
    pdb_file = relationship('PDBFile', primaryjoin="PDBMolecule.PDBFileID==PDBFile.ID")

    # Children relationships
    chains = relationship("PDBMoleculeChain", primaryjoin="and_(PDBMoleculeChain.PDBFileID==PDBMolecule.PDBFileID, PDBMoleculeChain.MoleculeID==PDBMolecule.MoleculeID)")

    def __repr__(self):
        return 'PDBMolecule ({0}-{1}). Name: {2}. Organism: {3}'.format(self.PDBFileID, self.MoleculeID, '/'.join([s for s in [self.Name or '', self.Synonym or ''] if s]), self.Organism or 'N/A')


class PDBMoleculeChain(DeclarativeBase):
    __tablename__ = 'PDBMoleculeChain'

    PDBFileID = Column(String(10), ForeignKey('PDBMolecule.PDBFileID'), ForeignKey('PDBChain.PDBFileID'), nullable=False, primary_key=True)
    MoleculeID = Column(Integer, ForeignKey('PDBMolecule.MoleculeID'), nullable=False, primary_key=True)
    Chain = Column(String(1), ForeignKey('PDBChain.Chain'), nullable=False, primary_key=True)

    # Parent relationships
    pdb_molecule = relationship('PDBMolecule', primaryjoin="and_(PDBMoleculeChain.PDBFileID==PDBMolecule.PDBFileID, PDBMoleculeChain.MoleculeID==PDBMolecule.MoleculeID)")
    pdb_chain = relationship('PDBChain', primaryjoin="and_(PDBMoleculeChain.PDBFileID==PDBChain.PDBFileID, PDBMoleculeChain.Chain==PDBChain.Chain)")

    def __repr__(self):
        return 'PDBMoleculeChain ({0}-{1}-{2}).'.format(self.PDBFileID, self.MoleculeID, self.Chain)


class PDBResidue(DeclarativeBase):
    __tablename__ = 'PDBResidue'

    ID = Column(Integer, nullable=False, primary_key=True)
    PDBFileID = Column(String(10), ForeignKey('PDBChain.PDBFileID'), nullable=False)
    Chain = Column(String(1), ForeignKey('PDBChain.Chain'), nullable=False)
    ResidueID = Column(String(5), nullable=False)
    ResidueAA = Column(String(1), ForeignKey('AminoAcid.Code'), nullable=False)
    ResidueType = Column(Enum('Protein', 'DNA', 'RNA'), nullable=False, default=u'Protein')
    IndexWithinChain = Column(Integer, nullable=False)
    CoordinatesExist = Column(TINYINT(1), nullable=False)
    RecognizedByRosetta = Column(TINYINT(1), nullable=True)
    BFactorMean = Column(DOUBLE, nullable=True)
    BFactorDeviation = Column(DOUBLE, nullable=True)
    MonomericExposure = Column(DOUBLE, nullable=True)
    MonomericDSSP = Column(String(1), nullable=True)
    ComplexExposure = Column(DOUBLE, nullable=True)
    ComplexDSSP = Column(String(1), nullable=True)

    # Parent relationships
    pdb_chain = relationship('PDBChain', primaryjoin="and_(PDBResidue.PDBFileID==PDBChain.PDBFileID, PDBResidue.Chain==PDBChain.Chain)")

    # Child relationships
    residue = relationship('AminoAcid', primaryjoin="PDBResidue.ResidueAA==AminoAcid.Code")

    def __repr__(self):
        return 'PDBResidue {0}. {1} {2} {3} ({4}). Exposure: {5}. DSSP: {6}.'.format(self.residue.LongCode, self.PDBFileID, self.Chain, (self.ResidueAA + self.ResidueID.strip()).ljust(5), self.ResidueType, self.ComplexExposure, self.ComplexDSSP)


class PDBLigand(DeclarativeBase):
    __tablename__ = 'PDBLigand'

    PDBFileID = Column(String(10), ForeignKey('PDBChain.PDBFileID'), nullable=False, primary_key=True)
    Chain = Column(String(1), ForeignKey('PDBChain.Chain'), nullable=False, primary_key=True)
    SeqID = Column(String(5), nullable=False, primary_key=True)
    PDBLigandCode = Column(String(3), nullable=False)
    LigandID = Column(Integer, ForeignKey('Ligand.ID'), nullable=False)
    ParamsFileContentID = Column(Integer, ForeignKey('FileContent.ID'), nullable=True)


class PDBIon(DeclarativeBase):
    __tablename__ = 'PDBIon'

    PDBFileID = Column(String(10), nullable=False, primary_key=True)
    Chain = Column(String(1), nullable=False, primary_key=True)
    SeqID = Column(String(5), nullable=False, primary_key=True)
    PDBIonCode = Column(String(3), nullable=False)
    IonID = Column(Integer, nullable=False)
    ParamsFileContentID = Column(Integer, nullable=True)
    Element = Column(String(2), nullable=False)


#########################################
#                                       #
#  Publications and associated records  #
#                                       #
#########################################


class Publication(DeclarativeBase):
    __tablename__ = 'Publication'

    ID = Column(String(64), nullable=False, primary_key=True)
    DGUnit = Column(Enum('kJ/mol','kcal/mol','cal/mol'), nullable=True)
    DDGConvention = Column(Enum('Rosetta','ProTherm','Unknown','Not applicable'), nullable=True)
    Notes = Column(Text, nullable=True)
    DGNotes = Column(Unicode(1024), nullable=True)
    DGUnitUsedInProTherm = Column(Enum('kcal/mol','kJ/mol'), nullable=True)
    DDGProThermSignNotes = Column(String(1024), nullable=True)
    DDGValuesNeedToBeChecked = Column(TINYINT(1), nullable=False, default=0)
    RIS = Column(Text, nullable=True)
    Title = Column(Unicode(256), nullable=True)
    Publication = Column(String(256), nullable=True)
    Volume = Column(String(8), nullable=True)
    Issue = Column(String(8), nullable=True)
    StartPage = Column(String(16), nullable=True)
    EndPage = Column(String(16), nullable=True)
    PublicationYear = Column(Integer, nullable=True)
    PublicationDate = Column(DateTime, nullable=True)
    DOI = Column(String(64), nullable=True)
    URL = Column(String(128), nullable=True)


class PublicationAuthor(DeclarativeBase):
    __tablename__ = 'PublicationAuthor'

    PublicationID = Column(String(64), nullable=False, primary_key=True)
    AuthorOrder = Column(Integer, nullable=False, primary_key=True)
    FirstName = Column(Unicode(64), nullable=False)
    MiddleNames = Column(Unicode(64), nullable=True)
    Surname = Column(Unicode(64), nullable=True)


class PublicationIdentifier(DeclarativeBase):
    __tablename__ = 'PublicationIdentifier'

    SourceID = Column(String(64), nullable=False, primary_key=True)
    ID = Column(String(256), nullable=False, primary_key=True)
    Type = Column(Enum('URL','DOI','ISSN','ESSN','PMID','MANUAL'), nullable=False)


class PublicationDDGValueLocation(DeclarativeBase):
    __tablename__ = 'PublicationDDGValueLocation'

    SourceID = Column(String(64), nullable=False, primary_key=True)
    Location = Column(String(256), nullable=False, primary_key=True)
    Notes = Column(String(512), nullable=True)


###########################
#                         #
#  Bookkeeping functions  #
#                         #
###########################




def generate_sqlalchemy_definition(tablenames = []):
    '''This function generates the SQLAlchemy class definitions from the database. The generation does not parse the
       entire definition - it omits unique keys, foreign key constraints etc. but it saves a lot of manual work setting
       up the boilerplate field definitions. When the database schema changes, call this function to update the
       SQLAlchemy class definitions. You may want/need to reuse any existing relationships defined between tables.'''
    sc = MySQLSchemaConverter('kortemmelab', 'kortemmelab.ucsf.edu', 'ddG', read_file(os.path.join('..', 'pw')).strip(), 3306, "/var/lib/mysql/mysql.sock")
    #sc.get_sqlalchemy_schema(['PDBFile', 'PDBChain', 'PDBMolecule', 'PDBMoleculeChain', 'PDBResidue'])
    sc.get_sqlalchemy_schema(tablenames)


def test_schema_against_database_instance(DDG_db):
    '''Make sure that our SQLAlchemy definitions match the database. This should be run by the API prior to connection
       as it lets the admin know that they need to update the schema here (use generate_sqlalchemy_definition to update
       the schema).'''
    database_to_class_mapping = {}
    clsmembers = inspect.getmembers(sys.modules[__name__], inspect.isclass)
    clsmembers = [c[1] for c in clsmembers if issubclass(c[1], DeclarativeBase) and c[1] != DeclarativeBase]
    for c in clsmembers:
        database_to_class_mapping[c.__tablename__] = c

    inconsistencies = []
    for tblname, pcls in sorted(database_to_class_mapping.iteritems()):
        represented_columns = set([c.name for c in list(sqlalchemy_inspect(pcls).columns)])
        tbl_columns = set([c['Field'] for c in DDG_db.execute_select('SHOW COLUMNS FROM {0}'.format(tblname))])
        if sorted(represented_columns) != sorted(tbl_columns):
            inconsistencies.append(pcls.__name__)
            colortext.error('The SQLAlchemy class {0} does not match the database schema.'.format(pcls.__name__))
            if represented_columns.difference(tbl_columns):
                colortext.warning('  The SQLAlchemy class has extra columns: {0}'.format(', '.join(sorted(represented_columns.difference(tbl_columns)))))
            if tbl_columns.difference(represented_columns):
                colortext.pcyan('  The MySQL schema definition has extra columns: {0}'.format(', '.join(sorted(tbl_columns.difference(represented_columns)))))
    if inconsistencies:
        generate_sqlalchemy_definition(inconsistencies)
        raise colortext.Exception('The following SQLAlchemy classes do not match the database schema: {0}.'.format(', '.join(inconsistencies)))



if __name__ == '__main__':
    generate_sqlalchemy_definition(['FileContent'])
    #generate_sqlalchemy_definition(['AminoAcid'])
    sys.exit(0)
    from ppi_api import get_interface as get_ppi_interface
    ppi_api = get_ppi_interface(read_file(os.path.join('..', 'pw')).strip())
    test_schema_against_database_instance(ppi_api.DDG_db)
