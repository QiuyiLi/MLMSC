import skbio
from .util import *


class TreeTableEntry:
    def __init__(self):
        self.__id = None
        self.__fakeId = None
        self.__name = None
        self.__parent = None
        self.__distanceToParent = None
        self.__children = []
        self.__distanceToChildren = []
        self.__clades = []
        self.__splits = []

    def __repr__(self):
        return (f'<TreeTableEntry, id: {self.__id}, fakeId: {self.__fakeId}, name: {self.__name}, '
                f'parent: {self.__parent}, distanceToParent: {self.__distanceToParent}, '
                f'children: {self.__children}, distanceToChildren: {self.__distanceToChildren}, '
                f'clades: {self.__clades}, splits: {self.__splits}>')

    def __str__(self):
        return (f'<TreeTableEntry, id: {self.__id}, fakeId: {self.__fakeId}, name: {self.__name}, '
                f'parent: {self.__parent}, distanceToParent: {self.__distanceToParent}, '
                f'children: {self.__children}, distanceToChildren: {self.__distanceToChildren}, '
                f'clades: {self.__clades}, splits: {self.__splits}>')

    @property
    def id(self):
        return self.__id
    @id.setter
    def id(self, id):
        self.__id = id

    @property
    def fakeId(self):
        return self.__fakeId
    @fakeId.setter
    def fakeId(self, fakeId):
        self.__fakeId = fakeId

    @property
    def name(self):
        return self.__name
    @name.setter
    def name(self, name):
        self.__name = name

    @property
    def parent(self):
        return self.__parent
    @parent.setter
    def parent(self, parent):
        self.__parent = parent

    @property
    def distanceToParent(self):
        return self.__distanceToParent
    @distanceToParent.setter
    def distanceToParent(self, distanceToParent):
        self.__distanceToParent = distanceToParent

    @property
    def children(self):
        return self.__children
    @children.setter
    def children(self, children):
        self.__children = children

    @property
    def distanceToChildren(self):
        return self.__distanceToChildren
    @distanceToChildren.setter
    def distanceToChildren(self, distanceToChildren):
        self.__distanceToChildren = distanceToChildren

    @property
    def clades(self):
        return self.__clades
    @clades.setter
    def clades(self, clades):
        self.__clades = clades

    @property
    def splits(self):
        return self.__splits
    @splits.setter
    def splits(self, splits):
        self.__splits = splits

class TreeTable:
    def __init__(self):
        self.__skbioTree = None
        self.__table = []
        self.__tableDictId = {}
        self.__tableDictName = {}
        self.__root = None
        self.__leaves = []
        self.__treeHeight = -1

    def __repr__(self):
        string = '<TreeTable, \n'
        for entry in self.__table:
            string += '  ' + str(entry) + '\n'
        string += '>'
        return string

    def __str__(self):
        string = '<TreeTable, \n'
        for entry in self.__table:
            string += '  ' + str(entry) + '\n'
        string += '>'
        return string

    @property
    def skbioTree(self):
        return self.__skbioTree
    @skbioTree.setter
    def skbioTree(self, skbioTree):
        leaveNode = None
        for node in skbioTree.tips():
            leaveNode = node
            break
        self.__treeHeight = skbioTree.distance(leaveNode)
        self.__skbioTree = skbioTree

    @property
    def table(self):
        return self.__table

    @property
    def root(self):
        return self.__root

    @property
    def leaves(self):
        return self.__leaves

    @property
    def treeHeight(self):
        treeHeight = self.__treeHeight
        if treeHeight < 0:
            treeHeight = 0
        return treeHeight

    def getEntryById(self, id):
        return self.__tableDictId[id]

    def getEntryByName(self, name):
        return self.__tableDictName[name]
    
    def getFakeIdFromId(self, id):
        return self.__tableDictId[id].fakeId

    def createFromEntries(self, entries, skbioTree):
        self.__skbioTree = skbioTree
        for entry in entries:
            self.__table.append(entry)
            self.__tableDictId[entry.id] = entry
            self.__tableDictName[entry.name] = entry
            if entry.name == skbioTree.name:
                self.__root = entry
            if skbioTree.find(entry.name).is_tip():
                self.__leaves.append(entry)
        
        # sort the table by id
        self.__table.sort(key=lambda x: x.id)

        # calculate tree height (any leaf node to the root)
        self.__treeHeight = self.__distanceToRoot(self.leaves[0].id)

    def createFromSkbioTree(self, skbioTree, rename=True):
        if rename:
            # rename all tree nodes
            self.__renameTreeNodes(skbioTree)

        # assign ids in reversed time order
        queue = Queue()
        visited = set()
        for treeNode in skbioTree.tips():
            queue.push(treeNode)
        i = 0
        while not queue.isEmpty():
            treeNode = queue.pop()
            visited.add(treeNode)
            treeNode.id = i
            i += 1
            if treeNode.is_root():
                continue  # equivalently break
            elif all(True if child in visited else False 
                     for child in treeNode.parent.children):
                queue.push(treeNode.parent)

        # create entry for each tree node and store in the table
        for treeNode in skbioTree.traverse():
            entry = TreeTableEntry()
            entry.id = treeNode.id
            entry.name = treeNode.name
            self.__tableDictId[entry.id] = entry
            self.__tableDictName[entry.name] = entry
            if not treeNode.parent:
                entry.parent = -1
                entry.distanceToParent = -1.0
            else:
                entry.parent = treeNode.parent.id
                entry.distanceToParent = treeNode.distance(treeNode.parent)
                parentTreeNode = self.__tableDictId[treeNode.parent.id]
                parentTreeNode.children.append(entry.id)
                parentTreeNode.distanceToChildren.append(
                    treeNode.distance(treeNode.parent))
            self.__table.append(entry)

            if treeNode.is_tip():
                self.__leaves.append(entry)

        # sort the table by id
        self.__table.sort(key=lambda x: x.id)

        # get the root
        self.__root = self.__table[-1]

        # assign fake ids in post order
        self.__assignFakeIds(skbioTree)

        # calculate tree height (any leaf node to the root)
        self.__treeHeight = self.__distanceToRoot(self.leaves[0].id)

        self.__skbioTree = skbioTree

        return skbioTree

    def createFromNewickFile(self, path):
        f = open(path)
        skbioTree = skbio.read(f, format="newick", into=skbio.tree.TreeNode)
        f.close()
        return self.createFromSkbioTree(skbioTree)

    def __renameTreeNodes(self, skbioTree):
        if skbioTree.name:
            return skbioTree.name
        else:
            name = ''
            for child in skbioTree.children:
                name += self.__renameTreeNodes(child)
            skbioTree.name = name
            return skbioTree.name

    def __assignFakeIds(self, skbioTree):
        index = 0
        for treeNode in skbioTree.postorder():
            entry = self.getEntryById(treeNode.id)
            entry.fakeId = index
            index += 1

    def __distanceToRoot(self, nodeId):
        """
        find the distance of a given node to the root 
        needed when finding the walking distance
        """
        if (self.getEntryById(nodeId).parent < 0 or 
            nodeId == self.root.id):
            return 0
        else:
            distanceToParent = self.getEntryById(nodeId).distanceToParent
            parent = self.getEntryById(nodeId).parent
            return distanceToParent + self.__distanceToRoot(parent)

    def distanceToLeaf(self, nodeId, branchDistance):
        """
        given a coalescent event happening at "branchDistance" above 
        a speices node with "nodeId" find the distance of this event 
        to the bottom of the tree needed when assigning ids to the 
        coalescent tree
        """
        distance = branchDistance + (self.treeHeight - self.__distanceToRoot(nodeId))
        if distance < 0:
            distance = 0

        return distance