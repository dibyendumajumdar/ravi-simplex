local luasimplex = require("luasimplex")
local rsm
if compiler and compiler.loadfile then
    rsm = compiler.loadfile('luasimplex/rsm.lua')()
else
    rsm = require("luasimplex.rsm")
end

local M =
{
  -- number of variables
  nvars = 4,
  -- number of constraints
  nrows = 2,
  indexes = @integer[] {1, 2, 3, 1, 2, 4},
  elements = @number[] {1, 2, 1, 2, 1, 1},
  row_starts = @integer[] {1, 4, 7},
  c = @number[] {-1, -1, 0, 0},
  xl = @number[] {0, 0, 0, 0},
  xu = @number[] {math.huge, math.huge, math.huge, math.huge},
  b = @number[] {3, 3},
}

local I = luasimplex.new_instance(M.nrows, M.nvars)
rsm.initialise(M, I, {})

objective, x = rsm.solve(M, I, {})

io.stderr:write(("Objective: %g\n"):format(objective))
io.stderr:write("  x:")
for i = 1, M.nvars do io.stderr:write((" %g"):format(x[i])) end
io.stderr:write("\n")


-- EOF -------------------------------------------------------------------------

