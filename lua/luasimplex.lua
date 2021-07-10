-- Utilities for the simplex method


local luasimplex = {}

-- RSM error -------------------------------------------------------------------

local rsm_exception = {}

function rsm_exception:new(error, M, I, S)
  local e = { error=error, M=M, I=I, S=S }
  self.__index = self
  return setmetatable(e, self)
end


function rsm_exception:__tostring()
  return self.error
end


function luasimplex.error(e, M, I, S)
  error(rsm_exception:new(e, M, I, S), 2)
end

function luasimplex.array_init(no_ffi)
end

function luasimplex.new_model(nrows: integer, nvars: integer, nonzeroes)
  local M
  M = {}

  M.nvars = nvars
  M.nrows = nrows
  M.nonzeroes = nonzeroes

  M.indexes = table.intarray(nonzeroes)
  M.row_starts = table.intarray(nrows+1)
  M.elements = table.numarray(nonzeroes)

  M.b = table.numarray(nrows)
  M.c = table.numarray(nvars)
  M.xl = table.numarray(nvars)
  M.xu = table.numarray(nvars)
  return M
end


function luasimplex.free_model(M)
end


function luasimplex.new_instance(nrows: integer, nvars: integer)
  local I = {}

  local total_vars: integer = nvars + nrows

  I.status = table.intarray(total_vars)
  I.basics = table.intarray(nrows)
  I.basic_cycles = table.intarray(nvars, 0)

  I.costs = table.numarray(nvars, 0)
  I.x = table.numarray(total_vars)
  I.xu = table.numarray(total_vars)
  I.xl = table.numarray(total_vars)

  I.basic_costs = table.numarray(nrows)
  I.pi = table.numarray(nrows, 0)
  I.reduced_costs = table.numarray(nvars, 0)
  I.gradient = table.numarray(nrows, 0)
  I.Binverse = table.numarray(nrows * nrows)

  I.TOLERANCE = 1e-7
  return I
end


function luasimplex.free_instance(I)
end


--------------------------------------------------------------------------------

return luasimplex


-- EOF -------------------------------------------------------------------------

