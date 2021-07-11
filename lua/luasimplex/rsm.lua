local math = require("math")
local error, type = error, type

local luasimplex = require("luasimplex")
local abs = math.abs

-- Constants -------------------------------------------------------------------

local TOLERANCE: number = 1e-7
local NONBASIC_LOWER: integer = 1
local NONBASIC_UPPER: integer = -1
local NONBASIC_FREE: integer = 2
local BASIC: integer = 0


-- Computation parts -----------------------------------------------------------

local function compute_pi(M: table, I: table)
  -- pi = basic_costs' * Binverse
  local nrows: integer, pi: number[], Bi: number[], TOL: number = M.nrows, I.pi, I.Binverse, I.TOLERANCE
  local basic_costs: number[] = I.basic_costs
  for i = 1, nrows do pi[i] = 0 end
  for i = 1, nrows do
    local c: number = basic_costs[i]
    if abs(c) > TOL then
      for j = 1, nrows do
        pi[j] = pi[j] + c * Bi[(i-1)*nrows + j]
      end
    end
  end
end


local function compute_reduced_cost(M: table, I: table)
  -- reduced cost = cost - pi' * A 
  local reduced_costs: number[], status: integer[], TOL: number = I.reduced_costs, I.status, I.TOLERANCE
  local indexes: integer[], elements: number[], row_starts: integer[] = M.indexes, M.elements, M.row_starts
  local pi: number[] = I.pi
  local nvars: integer = M.nvars
  local nrows: integer = M.nrows
  local costs: number[] = I.costs

  -- initialise with costs (phase 2) or zero (phase 1 and basic variables)
  for i = 1, nvars do
    reduced_costs[i] = status[i] ~= 0 and costs[i] or 0
  end

  -- Compute rcs 'sideways' - work through elements of A using each one once
  -- the downside is that we write to reduced_costs frequently
  for i = 1, nrows do
    local p = pi[i]
    if abs(p) > TOL then
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
  local TOL: number = -I.TOLERANCE
  -- Find the variable with the "lowest" reduced cost, keeping in mind that it might be at its upper bound

  local cycles: integer, minrc: number, entering_index: integer = math.maxinteger, 0.0, -1
  local nvars: integer = M.nvars
  local status: integer[] = I.status 
  local reduced_costs: number[] = I.reduced_costs
  local basic_cycles: integer[] = I.basic_cycles

  for i = 1, nvars do
    local s: integer, rc: number = status[i]
    if s == NONBASIC_FREE then
      rc = -abs(reduced_costs[i])
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
  local nrows: integer, Bi: number[] = M.nrows, I.Binverse
  local indexes: integer[], elements: number[], row_starts: integer[] = M.indexes, M.elements, M.row_starts
  local gradient: number[] = {}

  if g then
    gradient = g
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
  local TOL: number = I.TOLERANCE
  local status: integer[] = I.status
  local reduced_costs: number[] = I.reduced_costs
  local xu: number[] = I.xu
  local xl: number[] = I.xl
  local x: number[] = I.x
  local nrows: integer, nvars: integer = M.nrows, M.nvars
  local basics: integer[] = I.basics

  local s: integer = status[entering_index]
  if s == NONBASIC_FREE then
    s = reduced_costs[entering_index] > 0 and -1 or 1
  end

  local max_change: number, leaving_index: integer, to_lower = xu[entering_index] - xl[entering_index], -1

  for i = 1, nrows do
    local g: number = gradient[i] * -s
    if abs(g) > TOL then
      local j: integer, bound: number = basics[i]
      local found_bound: integer = 0

      if g > 0 then
        if xu[j] < math.huge then 
          bound = xu[j] 
          found_bound = 1
        end
      else
        if xl[j] > -math.huge then 
          bound = xl[j] 
          found_bound = 1
        end
      end

      if found_bound == 1 then
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
  local c: number = I.max_change
  local basics: integer[] = I.basics
  local x: number[] = I.x
  local gradient: number[] = I.gradient  
  local nrows: integer = M.nrows

  for i = 1, nrows do
    local j: integer = basics[i]
    x[j] = x[j] - c * gradient[i]
  end
end


local function update_Binverse(M: table, I: table)
  local nrows: integer, li: integer, Bi: number[] = M.nrows, I.leaving_index, I.Binverse
  local gradient: number[] = I.gradient  

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
  local nvars: integer = M.nvars
  local I_xu: number[] = I.xu
  local I_xl: number[] = I.xl
  local I_x: number[] = I.x
  local M_xu: number[] = M.xu
  local M_xl: number[] = M.xl
  local status: integer[] = I.status

  for ii = 1, nvars do
    local i: integer = ii + offset
    I_xu[i], I_xl[i] = M_xu[i], M_xl[i]
    if M_xl[i] == -math.huge and M_xu[i] == math.huge then
      I_x[i] = 0
      status[i] = NONBASIC_FREE
    elseif abs(M_xl[i]) < abs(M_xu[i]) then
      I_x[i] = M_xl[i]
      status[i] = NONBASIC_LOWER
    else
      I_x[i] = M_xu[i]
      status[i] = NONBASIC_UPPER
    end
  end
end


local function initialise_artificial_variables(M: table, I: table, offset: integer)
  local nrows: integer, nvars: integer = M.nrows, M.nvars
  local indexes: integer[], elements: number[], row_starts: integer[] = M.indexes, M.elements, M.row_starts
  local b: number[] = M.b
  local xu: number[] = I.xu
  local xl: number[] = I.xl
  local x: number[] = I.x
  local basic_costs: number[] = I.basic_costs
  local basics: integer[] = I.basics
  local status: integer[] = I.status

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

  local nrows: integer = M.nrows

  if not S.TOLERANCE then S.TOLERANCE = TOLERANCE end
  I.TOLERANCE = S.TOLERANCE

  initialise_real_variables(M, I, offset)
  initialise_artificial_variables(M, I, offset)

  local Binverse: number[] = I.Binverse
  for i = 1, nrows do Binverse[(i-1)*nrows + i + offset] = 1 end

  return I
end


-- Solve -----------------------------------------------------------------------

local function solve(M: table, I: table, S)
  local TOLERANCE: number = I.TOLERANCE
  local x: number[] = I.x
  local basic_costs: number[] = I.basic_costs
  local basics: integer[] = I.basics
  local status: integer[] = I.status
  local x: number[] = I.x
  local xu: number[] = I.xu
  local xl: number[] = I.xl
  local basic_cycles: integer[] = I.basic_cycles

  local nvars: integer, nrows: integer = M.nvars, M.nrows
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
          if basics[i] > nvars and abs(x[basics[i]]) > TOLERANCE  then
            luasimplex.error("Infeasible", M, I, S)
          end
        end
        I.costs = M.c
        for i = 1, nrows do
          if basics[i] <= nvars then
            basic_costs[i] = M.c[basics[i] ]
          end
        end
        I.phase = 2
      else
        break  -- optimal
      end
    else
      local entering_index: integer = I.entering_index

      basic_cycles[entering_index] = basic_cycles[entering_index] + 1

      compute_gradient(M, I, entering_index, I.gradient)
      local to_lower
      I.leaving_index, I.max_change, to_lower = find_leaving_variable(M, I, entering_index, I.gradient)
      if monitor then monitor(M, I, S, "leaving_variable") end

      local max_change: number = I.max_change
      local leaving_index: integer = I.leaving_index
      local costs: number[] = I.costs

      if I.phase == 2 and max_change >= math.huge / 2 then
        luasimplex.error("Unbounded", M, I, S)
      end

      if abs(max_change) > TOLERANCE then
        for i = 1, nvars do
          basic_cycles[i] = 0
        end
      end

      update_variables(M, I)
      x[entering_index] = x[entering_index] + max_change

      if leaving_index ~= -1 then
        update_Binverse(M, I)

        local rli: integer = basics[leaving_index]
        x[rli] = to_lower and xl[rli] or xu[rli]
        status[rli] = to_lower and NONBASIC_LOWER or NONBASIC_UPPER

        basics[leaving_index] = entering_index
        basic_costs[leaving_index] = costs[entering_index]

        status[entering_index] = BASIC
      else
        status[entering_index] = -status[entering_index]
      end
    end
  end

  local objective = 0
  for i = 1, nvars do
    objective = objective + x[i] * M.c[i]
  end

  return objective, x, I.iterations
end


--------------------------------------------------------------------------------
ravi.compile(initialise)
ravi.compile(initialise_artificial_variables)
ravi.compile(initialise_real_variables)
ravi.compile(compute_pi)
ravi.compile(compute_reduced_cost)
ravi.compile(compute_gradient)
ravi.compile(find_leaving_variable)
ravi.compile(find_entering_variable)
ravi.compile(update_Binverse)
ravi.compile(update_variables)
ravi.compile(solve)

return
{
  initialise            = initialise,
  solve                 = solve,
  compute_gradient      = compute_gradient,
  find_leaving_variable = find_leaving_variable,
}


-- EOF -------------------------------------------------------------------------

