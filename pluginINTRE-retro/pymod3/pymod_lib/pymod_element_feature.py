
class ElementFeature():
    def __init__(self, ID, name, start=1, end=1, description='', type_of_feature='None'):
        self.start = start
        self.end = end
        self.id = ID
        self.name = name
        self.type_of_feature = type_of_feature
        self.description = description
        self.parent_seq_start = None
        self.parent_seq_end = None

    def __repr__(self):
        return str(self.__dict__)

    def __deepcopy__(self, memodict={}):
        newfeature = self.__class__(self.id, self.name, self.start, self.end, self.description, self.type_of_feature)
        newfeature.__dict__.update(self.__dict__.copy())
        return newfeature
        # cls = self.__class__
        # result = cls.__new__(cls)
        # memodict[id(self)] = result
        # for k, v in self.__dict__.items():
        #     setattr(result, k, copy.deepcopy(v, memo))
        # return result

class DomainFeature(ElementFeature):
    def __init__(self, ID, name, start, end, evalue, color=None, description='', element=None):
        ElementFeature.__init__(self, ID, name, start, end, description, type_of_feature='domain')
        self.domain_color = color
        self.evalue = evalue
        self.element = element
        self.offset = None # symmetric, nterm offset and cterm offset if the seq is splitted into domains
        if self.element:
            self.sequence = self.element.my_sequence[int(start)-1:int(end)]
        else:
            self.sequence = None

    def set_element(self, element):
        self.element = element
        if self.element:
            self.sequence = self.element.my_sequence[int(self.start)-1:int(self.end)]
        else:
            self.sequence = None