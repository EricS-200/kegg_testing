from .atom import Atom

class Bond:
    def __init__(self, other_atom: Atom, type: str):
        self._other_atom = other_atom
        self._type = type

    @property
    def other_atom(self):
        return self._other_atom

    @property
    def type(self):
        return self._type

    @other_atom.setter
    def other_atom(self, other_atom):
        self._other_atom = other_atom

    @type.setter
    def type(self, type):
        self._type = type

    def __eq__(self, other):
        return self._other_atom == other.other_atom and self._type == other.type

    def __repr__(self):
        return f"Bond({self.other_atom}, {self.type})"