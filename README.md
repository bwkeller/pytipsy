pytipsy: low(ish) level tools to help out with TIPSY Files
==========================================================

pytipsy is a couple of handy helpers for working with GASOLINE2 and ChaNGa tipsy
outputs.  The module includes two functions and a special class:

* rtipsy: a rewrite of James Wadsley's IDL script for reading raw tipsy files
* wtipsy: the same as above, but for writing
* rarray: function for reading tipsy auxiliary arrays (for data not stored in
  the primary tipsy snapshot)
* warray: the same as above, but for writing
