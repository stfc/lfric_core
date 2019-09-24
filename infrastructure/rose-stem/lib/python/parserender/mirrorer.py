#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
'''
Upload a filetree to a server of some sort with arbitrary fileTransforms
en-route.
'''

from __future__ import print_function

from __future__ import absolute_import
from abc import ABCMeta, abstractmethod
import ftplib
import os
import os.path
import re
import netrc
import pwd
import shutil
import stat
import StringIO
import six.moves.urllib.parse
import xml.etree.ElementTree as ET
import six
from six import unichr

##############################################################################
class TreeVisitor(six.with_metaclass(ABCMeta)):
  @abstractmethod
  def visit( self, directory, subdirs, files ):
    '''
    Examines a directory.
    '''
    pass

##############################################################################
class GenerateTreeIndeces(TreeVisitor):
  def __init__( self, renderer ):
    self._renderer = renderer
    self._directoryList = {}
    self._indexFound = None

  def visit( self, directory, subdirs, files ):
    treeLevel = directory.count( os.sep )
    if treeLevel < self._indexFound:
      self._indexFound = None

    if 'index.html' in files:
      self._indexFound = treeLevel
      return

    if not self._indexFound:
      self._directoryList[directory] = subdirs + files

  def newFiles( self ):
    pageList = {}
    for filename, fileList in self._directoryList.items():
      content = self._renderer.render( '"{}" Directory'.format(filename), \
                                       fileList )
      pageList[os.path.join( filename, 'index.html' )] = content
    return pageList

##############################################################################
class XmlTransformation(six.with_metaclass(ABCMeta)):
  @abstractmethod
  def transform( self, root ):
    '''
    Transform the XML element tree in some way.
    '''
    pass

##############################################################################
class RemoveTagsWithClass(XmlTransformation):
  '''
  Removes any tag with the specified class from an XHTML document.
  '''
  def __init__( self, className ):
    self._className = className

  def transform( self, root ):
    # TODO: This does not handle the case where there is
    #       more than one class specified in the "class" attribute.
    xpath = ".//*[@class='{}']".format( self._className )
    for parent in root.findall( xpath + '/..' ):
      for element in parent.findall( xpath ):
        parent.remove( element )

##############################################################################
class AddRelativePrefix(XmlTransformation):
  '''
  Adds a prefix to all relative URLs found in HREF attributes.
  '''
  def __init__( self, prefix ):
    self._prefix = prefix
    if self._prefix[-1] != '/':
      self._prefix = self._prefix + '/'

  def transform( self, root ):
    xpath = './/*[@href]'
    for parent in root.findall( xpath + '/..' ):
      for element in parent.findall( xpath ):
        url = six.moves.urllib.parse.urlparse( element.attrib['href'] )
        if not url.scheme and not url.netloc:
          element.attrib['href'] = self._prefix + element.attrib['href']

##############################################################################
class NoNakedURLs(XmlTransformation):
  '''
  Any URLs which do not specify a file have a default index page added.
  This assumes that all files have an extension.
  '''
  def __init__( self, pageName ):
    self._defaultPageName = pageName

  def transform( self, root ):
    xpath = r'.//*[@href]'
    filePattern = re.compile( r'.*\..+$' )
    for parent in root.findall( xpath + '/..' ):
      for element in parent.findall( xpath ):
        if not filePattern.match( element.attrib['href'] ):
          if not element.attrib['href'].endswith( '/' ):
            element.attrib['href'] = element.attrib['href'] + '/'
          element.attrib['href'] = element.attrib['href'] + 'index.html'

##############################################################################
class Credentials(six.with_metaclass(ABCMeta)):
  @abstractmethod
  def getCredentials():
    '''
    Gets username and password.
    '''
    pass

##############################################################################
class NetrcCredentials(Credentials):
  def __init__( self, host ):
    '''
    Construct a Credentials object which will use a .netrc file as its source.

    host - Name of the host as it appears in .netrc
    '''
    self._host     = host

  def getCredentials( self ):
    '''
    Returns the first match found.

    The file is opened and read here rather than in the constructor to minimise
    the amount of time key matter spends in memory. Hopefully. This whole
    thing is a bit shoddy.
    '''
    credentials = netrc.netrc()
    login, account, password =  credentials.authenticators( self._host )
    return login, password

##############################################################################
class ObjectMissing(Exception):
  pass

##############################################################################
class Uploader(six.with_metaclass(ABCMeta)):
  @abstractmethod
  def upload( self, fileObject, filename ):
    '''
    Upload the contents of the fileObject to the server as filename.
    '''
    pass

  @abstractmethod
  def _removeDirectory( self, filename ):
    '''
    Removes the directory, if it exists.
    '''
    pass

  @abstractmethod
  def _renameDirectory( self, originalFilename, newFilename ):
    '''
    If the original filename exists, rename it to the new filename.
    '''
    pass

  @abstractmethod
  def _makeDirectory( self, filename ):
    '''
    Creates a directory on the mirror.
    '''
    pass

  def prepare( self, directory ):
    '''
    Prepares to upload to directory.
    '''
    self._workDirectory = directory + '.new'
    self._finalDirectory = directory
    self._previousDirectory = directory + '-previous'

    self._removeDirectory( self._workDirectory)
    self._makeDirectory( self._workDirectory )

  def ensureDirectory( self, filename ):
    '''
    Creates a directory on the mirror if none exists.
    '''
    self._makeDirectory( os.path.join( self._workDirectory, filename ) )

  def commit( self ):
    '''
    Replaces any existing file tree with the uploaded one.
    '''
    self._removeDirectory( self._previousDirectory )
    try: # self._finalDirectory may not exist
      self._renameDirectory( self._finalDirectory, self._previousDirectory )
    except ObjectMissing:
      pass
    self._renameDirectory( self._workDirectory, self._finalDirectory )

##############################################################################
class UploaderFile(Uploader):
  def __init__( self ):
    pass

  def upload( self, fileObject, filename ):
    absoluteFilename = os.path.join( self._workDirectory, filename )
    with open( absoluteFilename, 'wt' ) as destination:
      while True:
        block = fileObject.read( 1024 )
        if block == '':
          break
        destination.write( block )
    os.chmod( absoluteFilename, \
              stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IROTH )

  def _makeDirectory( self, filename ):
    if not os.path.exists( filename ):
      os.mkdir( filename )
      os.chmod( filename, \
                stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR \
                | stat.S_IRGRP | stat.S_IXGRP \
                | stat.S_IROTH | stat.S_IXOTH )

  def _removeDirectory( self, filename ):
    if os.path.exists( filename ):
      shutil.rmtree( filename )

  def _renameDirectory( self, originalFilename, newFilename ):
    if os.path.exists( originalFilename ):
      os.rename( originalFilename, newFilename )

##############################################################################
class UploaderFTP(Uploader):
  def __init__( self, host, credentials=None ):
    self._client = ftplib.FTP( host, timeout=30 )

    if credentials:
      username, password = credentials.getCredentials()
      self._client.login( username, password )
      del password

  def upload( self, fileObject, filename ):
    try:
      self._client.storbinary( 'STOR {}/{}'.format( self._workDirectory, \
                                                    filename ), fileObject )
    except Exception as ex:
      message = 'Unable to upload file to server: {}'
      raise Exception( message.format( filename ) )

  def _makeDirectory( self, filename ):
    objects = self._client.nlst( filename )
    if objects == []:
      self._client.mkd( filename )

  def _removeDirectory( self, filename ):
    objects = self._client.nlst( filename )
    if objects == []:
      # Object may not exist or may be an empty directory
      try:
        self._client.rmd( filename )
      except Exception:
        pass
      return

    if objects == [filename]:
      try:
        self._client.delete( filename )
      except Exception as ex:
        message = 'Unable to delete file "{}" on server: {}'
        raise Exception( message.format( filename, ex ) )
    else:
      for fobject in objects:
        self._removeDirectory( fobject )
      try:
        self._client.rmd( filename )
      except Exception as ex:
        message = 'Unable to delete directory "{}" on server: {}'
        raise Exception( message.format( filename, ex ) )

  def _renameDirectory( self, originalFilename, newFilename ):
    objects = self._client.nlst( originalFilename )
    if objects == []:
      raise ObjectMissing( 'No such file: {}'.format( originalFilename) )

    try:
      self._client.rename( originalFilename, newFilename )
    except Exception as ex:
      message = 'Unable to rename "{}" as "{}": {}'
      raise Exception( message.format( originalFilename, newFilename, ex ) )

##############################################################################
class Mirrorer:
  '''
  Mirror a web site to a remote server.
  '''
  def __init__( self, destination, fileTransforms=[], treeTransforms=[] ):
    '''
    Constructor taking upload URL and fileTransforms to apply to HTML files.
    '''
    url = six.moves.urllib.parse.urlparse( destination )

    ET.register_namespace( '', 'http://www.w3.org/1999/xhtml' )

    if url.scheme == 'file':
      self._uploader = UploaderFile()
      self._uploader.prepare( url.path )
    elif url.scheme == 'ftp':
      # Slice leading slash off path to make it relative.
      self._uploader = UploaderFTP( url.netloc, \
                                    NetrcCredentials( url.netloc ) )
      self._uploader.prepare( url.path[1:] )
    else:
      raise NotImplementedError( 'Only file and FTP destinations are supported' )

    self._fileTransforms = fileTransforms
    self._treeTransforms = treeTransforms

  def mirror( self, source ):
    ET.register_namespace( '', 'http://www.w3.org/1999/xhtml' )
    for root, dirs, files in os.walk( source ):
      for transformation in self._treeTransforms:
        transformation.visit( os.path.relpath( root, source ), dirs, files )

      for directory in dirs:
        absoluteFilename = os.path.join( root, directory )
        relativeFilename = os.path.relpath( absoluteFilename, source )
        self._uploader.ensureDirectory( relativeFilename )

      for filename in files:
        absoluteFilename = os.path.join( root, filename )
        with open( absoluteFilename, 'rt' ) as fileStream:
          transformedStream = self._transformFile( absoluteFilename, \
                                                  fileStream )
          self._uploader.upload( transformedStream, \
                                os.path.relpath( absoluteFilename, source ) )
          transformedStream.close()

    for transformation in self._treeTransforms:
      for filename, content in transformation.newFiles().items():
        fakeFile = StringIO.StringIO()
        print( content, file=fakeFile )
        fakeFile.seek( 0, os.SEEK_SET )
        transformedStream = self._transformFile( filename , fakeFile )
        self._uploader.upload( transformedStream, filename )
        transformedStream.close()
        fakeFile.close()

    self._uploader.commit()

  def _transformFile( self, filename, stream ):
    base, extension = os.path.splitext( filename )
    if extension in ['.html', '.xhtml', '.htm']:
      xmlParser = ET.XMLParser()
      # In the absense of a comprehensive list to just lift up and use we
      # include only those actually appearing in our documents.
      xmlParser.entity['nbsp']  = unichr(160)
      xmlParser.entity['sect']  = unichr(167)
      xmlParser.entity['ndash'] = unichr(2013)
      xmlParser.entity['mdash'] = unichr(2014)
      try:
        tree = ET.parse( stream, xmlParser )
      except ET.ParseError as ex:
        message = 'Failed to parse {}: {}'
        raise Exception( message.format( filename, ex ) )
      treeRoot = tree.getroot()

      for transformation in self._fileTransforms:
        transformation.transform( treeRoot )

      content = StringIO.StringIO()
      tree.write( content, method='html' )
      content.seek( 0, os.SEEK_SET )

      return content
    else:
      return stream
