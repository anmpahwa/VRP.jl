# Instance files

Each instance folder contains the following files

- arcs.csv: matrix of arc lengths

- customer_nodes.csv: list of pcikup and delivery customer nodes in the network
    - in: customer node index
    - jn: associated pickup or delivery customer node
    - x: location on the x-axis
    - y: location on the y-axis
    - qc: demand
    - tc: service time (duration)
    - te: earliest service time
    - tl: latest service time

- depot_nodes.csv: list of depot nodes in the network
    - in: depot node index
    - x: location on the x-axis
    - y: location on the y-axis
    - ts: working-hours start time
    - te: working-hours end time
    - co: operational cost per package
    - cf: fixed cost

- fuelstation_nodes.csv: list of fuel station nodes in the network
    - in: fuel station node index
    - jn: associated vehicle fuel-type
    - x: location on the x-axis
    - y: location on the y-axis
    - tv: re-fueling rate
    - co: operational cost per unit fuel consumed re-fueling
    - cf: fixed cost

- vehicles.csv: list of vehicles at the depot nodes
    - iv: vehicle index
    - jv: vehicle type index
    - jf: vehicle fuel type index
    - id: depot node index
    - qv: capacity
    - lv: range
    - sv: speed
    - tf: re-fueling time at the depot node
    - td: service time per package at the depot node
    - tc: parking time at a customer node
    - tw: driver working-hours
    - r: driver work-load (maximum vehicle-routes)
    - cd: distance-based operational cost
    - ct: time-based operational cost
    - cf: fixed cost

