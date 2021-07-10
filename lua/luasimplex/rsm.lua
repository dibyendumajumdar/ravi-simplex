local math = require("math")
local error, type = error, type

local luasimplex = require("luasimplex")


-- Constants -------------------------------------------------------------------

local TOLERANCE: number = 1e-7
local NONBASIC_LOWER: integer = 1
local NONBASIC_UPPER: integer = -1
local NONBASIC_FREE: integer = 2
local BASIC: integer = 0


-- Computation parts -----------------------------------------------------------

local function compute_pi(M: table, I: table)
  -- pi = basic_costs' * Binverse
  local nrows: integer, pi: number[], Bi: number[], TOL: number = @integer( M.nrows ), @number[]( I.pi ), @number[]( I.Binverse ), @number( I.TOLERANCE )
  local basic_costs: number[] = @number[]( I.basic_costs )
  for i = 1, nrows do pi[i] = 0 end
  for i = 1, nrows do
    local c: number = basic_costs[i]
    if math.abs(c) > TOL then
      for j = 1, nrows do
        pi[j] = pi[j] + c * Bi[(i-1)*nrows + j]
      end
    end
  end
end


local function compute_reduced_cost(M: table, I: table)
  -- reduced cost = cost - pi' * A 
  local reduced_costs: number[], status: integer[], TOL: number = @number[]( I.reduced_costs ), @integer[]( I.status ), @number( I.TOLERANCE )
  local indexes: integer[], elements: number[], row_starts: integer[] = @integer[]( M.indexes ), @number[]( M.elements ), @integer[]( M.row_starts )
  local pi: number[] = @number[]( I.pi )
  local nvars: integer = @integer( M.nvars )
  local nrows: integer = @integer( M.nrows )
  local costs: number[] = @number[]( I.costs )

  -- initialise with costs (phase 2) or zero (phase 1 and basic variables)
  for i = 1, nvars do
    reduced_costs[i] = status[i] ~= 0 and costs[i] or 0
  end

  -- Compute rcs 'sideways' - work through elements of A using each one once
  -- the downside is that we write to reduced_costs frequently
  for i = 1, nrows do
    local p = pi[i]
    if math.abs(p) > TOL then
      for j = row_starts[i], row_starts[i+1]-1 do
        local k: integer = indexes[j]
        if status[k] ~= 0 then
          reduced_costs[k] = reduced_costs[k] - p * elements[j]
        end        
      end
    end
  end
end


local function find_entering_variable(M: table, I: table)
  local TOL: number = -@number(I.TOLERANCE)
  -- Find the variable with the "lowest" reduced cost, keeping in mind that it might be at its upper bound

  local cycles: integer, minrc: number, entering_index: integer = @integer( math.maxinteger ), 0.0, -1
  local nvars: integer = @integer( M.nvars )
  local status: integer[] = @integer[]( I.status ) 
  local reduced_costs: number[] = @number[]( I.reduced_costs )
  local basic_cycles: integer[] = @integer[]( I.basic_cycles )

  for i = 1, nvars do
    local s: integer, rc: number = status[i]
    if s == NONBASIC_FREE then
      rc = - @number( math.abs(reduced_costs[i]) )
    else
      rc = s * reduced_costs[i]
    end
    local c: integer = basic_cycles[i]
    if (c < cycles and rc < TOL) or (c == cycles and rc < minrc) then
      minrc = rc
      cycles = basic_cycles[i]
      entering_index = i
    end
  end
  return entering_index
end


local function compute_gradient(M: table, I: table, entering_index: integer, g)
  -- gradient = Binverse * entering column of A
  local nrows: integer, Bi: number[] = @integer( M.nrows ), @number[]( I.Binverse )
  local indexes: integer[], elements: number[], row_starts: integer[] = @integer[]( M.indexes ), @number[]( M.elements ), @integer[]( M.row_starts )
  local gradient: number[] = {}

  if g then
    gradient = @number[]( g )
    for i = 1, nrows do gradient[i] = 0 end
  else
    gradient = table.numarray(nrows, 0)
  end

  for i = 1, nrows do
    local v: number
    local found: integer = 0
    for j = row_starts[i], row_starts[i+1]-1 do
      local column: integer = indexes[j]
      if column == entering_index then
        v = elements[j]
        found = 1
        break
      elseif column > entering_index then
        break
      end
    end
    if found == 1 then
      for j = 1, nrows do
        gradient[j] = gradient[j] + v * Bi[(j-1)*nrows + i]
      end
    end
  end
  return gradient
end


local function find_leaving_variable(M: table, I: table, entering_index: integer, gradient: number[])
  local TOL: number = @number( I.TOLERANCE )
  local status: integer[] = @integer[]( I.status )
  local reduced_costs: number[] = @number[]( I.reduced_costs )
  local xu: number[] = @number[]( I.xu )
  local xl: number[] = @number[]( I.xl )
  local x: number[] = @number[]( I.x )
  local nrows: integer, nvars: integer = @integer( M.nrows ), @integer( M.nvars )
  local basics: integer[] = @integer[]( I.basics )

  local s: integer = status[entering_index]
  if s == NONBASIC_FREE then
    s = reduced_costs[entering_index] > 0 and -1 or 1
  end

  local max_change: number, leaving_index: integer, to_lower = xu[entering_index] - xl[entering_index], -1

  for i = 1, nrows do
    local g: number = gradient[i] * -s
    if math.abs(g) > TOL then
      local j: integer, bound: number = basics[i]
      local found_bound = false

      if g > 0 then
        if xu[j] < math.huge then 
          bound = xu[j] 
          found_bound = true
        end
      else
        if xl[j] > -math.huge then 
          bound = xl[j] 
          found_bound = true
        end
      end

      if found_bound then
        local z: number = (bound - x[j]) / g
        -- we prefer to get rid of artificials when we can
        if z < max_change or (j > nvars and z <= max_change) then
          max_change = z
          leaving_index = i
          to_lower = g < 0
        end
      end
    end
  end
  
  return leaving_index, max_change * s, to_lower
end


local function update_variables(M: table, I: table)
  local c: number = @number( I.max_change )
  local basics: integer[] = @integer[]( I.basics )
  local x: number[] = @number[]( I.x )
  local gradient: number[] = @number[]( I.gradient )  
  local nrows: integer = @integer( M.nrows )

  for i = 1, nrows do
    local j: integer = basics[i]
    x[j] = x[j] - c * gradient[i]
  end
end


local function update_Binverse(M: table, I: table)
  local nrows: integer, li: integer, Bi: number[] = @integer( M.nrows ), @integer( I.leaving_index ), @number[]( I.Binverse )
  local gradient: number[] = @number[]( I.gradient )  

  local ilg: number = 1.0 / gradient[li]
  for i = 1, nrows do
    if i ~= li then
      local gr: number = gradient[i] * ilg
      for j = 1, nrows do
        Bi[(i-1)*nrows + j] = Bi[(i-1)*nrows + j] - gr * Bi[(li-1)*nrows + j]
      end
    end
  end
  for j = 1, nrows do
    Bi[(li-1)*nrows + j] = Bi[(li-1)*nrows + j] * ilg
  end
end


-- Initialisation --------------------------------------------------------------

local function initialise_real_variables(M: table, I: table, offset: integer)
  local nvars: integer = @integer( M.nvars )
  local I_xu: number[] = @number[]( I.xu )
  local I_xl: number[] = @number[]( I.xl )
  local I_x: number[] = @number[]( I.x )
  local M_xu: number[] = @number[]( M.xu )
  local M_xl: number[] = @number[]( M.xl )
  local status: integer[] = @integer[]( I.status )

  for ii = 1, nvars do
    local i: integer = ii + offset
    I_xu[i], I_xl[i] = M_xu[i], M_xl[i]
    if M_xl[i] == -math.huge and M_xu[i] == math.huge then
      I_x[i] = 0
      status[i] = NONBASIC_FREE
    elseif math.abs(M_xl[i]) < math.abs(M_xu[i]) then
      I_x[i] = M_xl[i]
      status[i] = NONBASIC_LOWER
    else
      I_x[i] = M_xu[i]
      status[i] = NONBASIC_UPPER
    end
  end
end


local function initialise_artificial_variables(M: table, I: table, offset: integer)
  local nrows: integer, nvars: integer = @integer( M.nrows ), @integer( M.nvars )
  local indexes: integer[], elements: number[], row_starts: integer[] = @integer[]( M.indexes ), @number[]( M.elements ), @integer[]( M.row_starts )
  local b: number[] = @number[]( M.b )
  local xu: number[] = @number[]( I.xu )
  local xl: number[] = @number[]( I.xl )
  local x: number[] = @number[]( I.x )
  local basic_costs: number[] = @number[]( I.basic_costs )
  local basics: integer[] = @integer[]( I.basics )
  local status: integer[] = @integer[]( I.status )

  for ii = 1, nrows do
    local i: integer = ii + offset
    local z: number = b[i]
    for j = row_starts[i], row_starts[i+1]-1 do
      z = z - elements[j] * x[indexes[j]]
    end
    local k: integer = nvars + i
    x[k] = z
    status[k] = BASIC
    basics[i] = k
    if z < 0 then
      basic_costs[i], xl[k], xu[k] = -1, -math.huge, 0
    else
      basic_costs[i], xl[k], xu[k] = 1, 0, math.huge
    end
    if type(M) == "table" and M.variable_names and M.constraint_names then
      M.variable_names[k] = M.constraint_names[i].."_ARTIFICIAL"
    end
  end
end


local function initialise(M: table, I: table, S, c_arrays)
  local offset: integer = c_arrays and -1 or 0

  local nrows: integer = @integer( M.nrows )

  if not S.TOLERANCE then S.TOLERANCE = TOLERANCE end
  I.TOLERANCE = S.TOLERANCE

  initialise_real_variables(M, I, offset)
  initialise_artificial_variables(M, I, offset)

  local Binverse: number[] = @number[]( I.Binverse )
  for i = 1, nrows do Binverse[(i-1)*nrows + i + offset] = 1 end

  return I
end


-- Solve -----------------------------------------------------------------------

local function solve(M, I, S)
  local TOLERANCE = I.TOLERANCE
  
  local nvars, nrows = M.nvars, M.nrows
  I.iterations = 0
  I.phase = 1
  local monitor = S.monitor

  while true do
    I.iterations = I.iterations + 1
    if monitor then monitor(M, I, S, "iteration") end
    if I.iterations > 10000 then
      luasimplex.error("Iteration limit", M, I, S)
    end

    compute_pi(M, I)
    compute_reduced_cost(M, I)
    I.entering_index = find_entering_variable(M, I)
    if monitor then monitor(M, I, S, "entering_variable") end

    if I.entering_index == -1 then
      if I.phase == 1 then
        for i = 1, nrows do
          if I.basics[i] > nvars and math.abs(I.x[I.basics[i] ]) > TOLERANCE  then
            luasimplex.error("Infeasible", M, I, S)
          end
        end
        I.costs = M.c
        for i = 1, nrows do
          if I.basics[i] <= nvars then
            I.basic_costs[i] = M.c[I.basics[i] ]
          end
        end
        I.phase = 2
      else
        break  -- optimal
      end
    else

      I.basic_cycles[I.entering_index] = I.basic_cycles[I.entering_index] + 1

      compute_gradient(M, I, I.entering_index, I.gradient)
      local to_lower
      I.leaving_index, I.max_change, to_lower = find_leaving_variable(M, I, I.entering_index, I.gradient)
      if monitor then monitor(M, I, S, "leaving_variable") end

      if I.phase == 2 and I.max_change >= math.huge / 2 then
        luasimplex.error("Unbounded", M, I, S)
      end

      if math.abs(I.max_change) > TOLERANCE then
        for i = 1, nvars do
          I.basic_cycles[i] = 0
        end
      end

      update_variables(M, I)
      I.x[I.entering_index] = I.x[I.entering_index] + I.max_change

      if I.leaving_index ~= -1 then
        update_Binverse(M, I)

        local rli = I.basics[I.leaving_index]
        I.x[rli] = to_lower and I.xl[rli] or I.xu[rli]
        I.status[rli] = to_lower and NONBASIC_LOWER or NONBASIC_UPPER

        I.basics[I.leaving_index] = I.entering_index
        I.basic_costs[I.leaving_index] = I.costs[I.entering_index]

        I.status[I.entering_index] = BASIC
      else
        I.status[I.entering_index] = -I.status[I.entering_index]
      end
    end
  end

  local objective = 0
  for i = 1, nvars do
    objective = objective + I.x[i] * M.c[i]
  end

  return objective, I.x, I.iterations
end


--------------------------------------------------------------------------------

return
{
  initialise            = initialise,
  solve                 = solve,
  compute_gradient      = compute_gradient,
  find_leaving_variable = find_leaving_variable,
}


-- EOF -------------------------------------------------------------------------

