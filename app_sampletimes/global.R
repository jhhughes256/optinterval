# set up times to be sourced by selectizeInput
  glob.times <- as.character(
    sort(
      unique(
        c(seq(from = 0, to = 5, by = 0.25), seq(from = 5, to = 24, by = 1))
      )
    )
  )

# set up default times
  def.times <- c(0, 0.25, 0.5, 1, 2, 4, 8, 12, 24)
