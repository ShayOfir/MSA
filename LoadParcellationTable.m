function ParcTab = LoadParcellationTable (fn, sheet)

ParcTab = readtable(fn, 'Sheet', sheet);

end