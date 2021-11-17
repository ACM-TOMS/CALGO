package Symbol_Table is

-- DESCRIPTION :
--   This package supplies a data-abstraction to work with a symbol table.

-- AUXILIARY DATA STRUCTURES :

  subtype symbol is string(1..3);

-- EXCEPTIONS :
    
    OVERFLOW_IN_THE_SYMBOL_TABLE : exception;
     -- occurs when a new symbol is added to a full symbol table

    INDEX_OUT_OF_RANGE : exception;
     -- occurs when a symbol is asked that is not in the range of the table

-- CREATORS :

  procedure Init ( max : in natural );

  -- DESCRIPTION :
  --   A new symbol table is created with place for max symbols.

  procedure Enlarge ( max : in natural );

  -- DESCRIPTION :
  --   Enlarges the symbol table so that it can contain as many symbols
  --   as the previous maximum plus the new max.

-- CONSTRUCTOR :

  procedure Add ( sb : in symbol );
  procedure Add ( sb : in symbol; pos : out natural );

  -- DESCRIPTION :
  --   A new symbol is added to the symbol table;
  --   pos is the entry of the added symbol in the table.

-- SELECTORS :

  function "<" ( s1,s2 : Symbol ) return boolean;
  function ">" ( s1,s2 : Symbol ) return boolean;

  -- DESCRIPTION :
  --   Defines an order relation on the symbols.

  function Number return natural;
    
  -- DESCRIPTION :
  --   Returns the number of current symbols in the table.

  function Empty return boolean;

  -- DESCRIPTION :
  --   Returns true if the symbol table has not been initialized yet,
  --   or if a Clear has been done.

  function Get ( sb : symbol ) return natural;

  -- DESCRIPTION :
  --  The entry of the symbol in the table is returned.
  --  If the symbol does not occur in the table, then 0 is returned.

  function Get ( i : natural ) return symbol;

  -- DESCRIPTION :
  --   The symbol corresponding with the ith unknown is returned.

-- DESTRUCTOR :

  procedure Clear;
  
  -- DESCRIPTION :
  --   The allocated memory space is freed.

end Symbol_Table;
