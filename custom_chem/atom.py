class Atom:
    def __init__(self, mapping_num: int, symbol: str, kegg_atom_type: str = None, charge: int = 0):
        self._mapping_num = mapping_num
        self._symbol = symbol
        self._kegg_atom_type = kegg_atom_type
        self._charge = charge

    @property
    def mapping_num(self):
        return self._mapping_num

    @property
    def symbol(self):
        return self._symbol

    @property
    def kegg_atom_type(self):
        return self._kegg_atom_type

    @property
    def charge(self):
        return self._charge

    @mapping_num.setter
    def mapping_num(self, mapping_num):
        self._mapping_num = mapping_num

    @symbol.setter
    def symbol(self, symbol):
        self._symbol = symbol

    @kegg_atom_type.setter
    def kegg_atom_type(self, kegg_atom_type):
        self._kegg_atom_type = kegg_atom_type
        self._charge = 0

    @charge.setter
    def charge(self, charge):
        self._charge = charge

    def __eq__(self, other):
        return self._mapping_num == other.mapping_num and self._symbol == other.symbol

    def __hash__(self):
        return hash((self._mapping_num, self._symbol))

    def __repr__(self):
        return f"Atom({self._mapping_num}, {self._symbol}, {self._kegg_atom_type}, {self._charge})"