# TODO
- Reduce Solution datastructure (remove φ)
- Modify insertnode!, removenode!, movevehicle!, worstcustomer!, intermove!, intraswap!, interswap!, and interopt!
- Add refueling stations (F = {f; f ∈ F})
- Update customer node datastructure to include the nearest refueling station (f)
- Update build, isfeasible, insertnode!, removenode!, and movevehicle! (do update customer node indices if a refueling station is visited)
- When updating en-route variables for a particular route, all visits to the refueling stations must be removed for this and subsequent routes, and then re-inserted appropriately