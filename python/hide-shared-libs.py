#!/usr/bin/env python

# Author: Todd Kulesza <todd@dropline.net>
# Released to the public domain

import argparse
import os
import shutil

class Hider:
    args = None
    rootDir = None
    srcSuffix = None
    dstSuffix = None

    def run(self):
        self.args = self.parseArgs()
        if self.args.hide:
            self.srcSuffix = ".dylib"
            self.dstSuffix = ".hidden"
        else:
            self.srcSuffix = ".hidden"
            self.dstSuffix = ".dylib"
        libs = self.listSourceFiles()
        self.moveFiles(libs)

    def parseArgs(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('-d', '--dir', help="Root directory", required=True)
        parser.add_argument('-m', '--hide', help="Hide shared libraries", action='store_true')
        parser.add_argument('-r', '--restore', help="Restore shared libraries", action='store_true')
        return parser.parse_args()

    def listSourceFiles(self):
        libraries = []
        for path, subdirs, files in os.walk(self.args.dir):
            for name in files:
                if self.isSourceFile(name):
                    libraries.append(os.path.join(path, name))
        libraries.sort()
        return libraries

    def isSourceFile(self, library):
        return library.endswith(self.srcSuffix)

    def moveFiles(self, srcFiles):
        for path in srcFiles:
            dstPath = path.replace(self.srcSuffix, self.dstSuffix)
            print("Moving '{}' to '{}'".format(path, dstPath))
            shutil.move(path, dstPath)

if __name__ == "__main__":
    hider = Hider()
    hider.run()
