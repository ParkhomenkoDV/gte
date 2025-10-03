## gte.nodes

gte.nodes = дирeктория с узлами ГТД

- Compressor

```
          ------------
         |            |
inlet -> | Compressor | -> outlet
         |            |
          ------------
               
```

- CombustionChamber

```
                 fuel
                   |
                   v
          -------------------
         |                   |
inlet -> | CombustionChamber | -> outlet
         |                   |
          -------------------
```

- Turbine

```
           cooling
              |
              v
          ---------
         |         |
inlet -> | Turbine | -> outlet
         |         |
          ---------
```

## Shaft

```
Compressor_1 ... Compressor_n Turbine_1 ... Turbine_n Load_1 ... Load_n
     |                |           |             |        |          |
      --------------------------------------------------------------
                                   Shaft
```

## Mixing

```
            ----------
inlet_1 -> |          |
inlet_2 -> |  Mixing  | -> outlet
inlet_n -> |          |
            ----------
```