#!/usr/bin/awk -f

# postdiff.awk -- postprocessor for checking outputs of PWscf examples
# checking is done in three steps: preprocess, diff against reference
# data, postprocess
# this way "false positives" are eliminated, that is, differences that
# don't mean that something went wrong

# for each (group of) line(s)
{
  # read block of lines from output of diff
  read_diff(); # this sets s1, s2, n1, n2, line1, line2

  # check whether block can be ignored, if not, print it
  while (n1 > 0 || n2 > 0)
    check_diff();
}

function read_diff()
{
  # read output of diff
  read_head(); # this sets s1, s2, n1, n2
  for (i=0; i<n1; i++)
    {
      getline; line1[i] = $0;
    }
  if (n1 > 0 && n2 > 0)
    getline; # separator
  for (i=0; i<n2; i++)
    {
      getline; line2[i] = $0;
    }
}

function check_diff()
{
  i1 = i2 = 0;

  # extract block of lines with the same key
  if (n1 > 0)
    {
      key1 = get_key(line1[0]);
      for (i1=1; i1<n1 && get_key(line1[i1]) == key1; i1++)
	;
    }
  if (n2 > 0)
    {
      key2 = get_key(line2[0]);
      for (i2=1; i2<n2 && get_key(line2[i2]) == key2; i2++)
	;
    }

  # choose which key comes first
  if (i1 == 0)
    key = key2;
  else if (i2 == 0)
    key = key1;
  else if (key1 != key2)
    {
      # look for matching keys
      for (j1=i1; j1<n1 && get_key(line1[j1]) != key2; j1++)
	;
      for (j2=i2; j2<n2 && get_key(line2[j2]) != key1; j2++)
	;

      if (j1 == n1 && j2 < n2)
	{
	  # match found for key1, do key2 first
	  key = key2; i1 = 0;
	}
      else if (j2 == n2 && j1 < n1)
	{
	  # match found for key2, do key1 first
	  key = key1; i2 = 0;
	}
      else
	{
	  # either no key or both keys match
	  # heuristics: larger block of lines comes first
	  if (n1 >= n2)
	    {
	      key = key1; i2 = 0;
	    }
	  else
	    {
	      key = key2; i1 = 0;
	    }
	}
    }
  else
    key = key1;

  # recreate output of diff
  if (! check_key())
    print_diff();

  # cut away the block of lines just examined
  s1 += i1;
  s2 += i2;
  n1 -= i1;
  n2 -= i2;
  for (k=0; k<n1; k++)
    line1[k] = line1[k + i1];
  for (k=0; k<n2; k++)
    line2[k] = line2[k + i2];
}

function check_key()
{
  if (key == "TIMING")
    {
      return 1;
    }
  else if (key == "PARALLEL")
    {
      # always accept
      return 1;
    }
  else if (key == "FFTGRID")
    {
      # number of lines must match
      if (i1 != i2)
	return 0;

      # check passed
      return 1;
    }
  else if (key == "KWEIGHT")
    {
      tolerance = 2e-7;

      # number of lines must match
      if (i1 != i2)
	return 0;

      # all pairs of lines must match
      for (j=0; j<i1; j++)
	{
	  split(line1[j], x1);
	  split(line2[j], x2);
	  delta = x1[12] - x2[12];
	  if (delta*delta > tolerance*tolerance)
	    return 0;
	}

      # check passed
      return 1;
    }
  else if (key == "ITERATION")
    {
      # discard everything
      return 1;
    }
  else if (key == "FERMI")
    {
      tolerance = 5e-4;

      # there must be exactly one line
      if (i1 != 1 || i2 != 1)
	return 0;

      split(line1[0], x1);
      split(line2[0], x2);
      delta = x1[7] - x2[7];
      if (delta*delta > tolerance*tolerance)
	return 0;

      # check passed
      return 1;
    }
  else if (key == "ENERGY")
    {
      tolerance = 1e-5;

      # number of lines must match
      if (i1 != i2)
	return 0;

      # all pairs of lines must match
      for (j=0; j<i1; j++)
	{
	  # energy is the field after = or < (not the < prepended by diff!)
	  split(line1[j], x1);
	  split(line2[j], x2);
	  for (k in x1)
	    if (k > 2 && match(x1[k], "[=<]"))
	      break;
	  delta = x1[k+1] - x2[k+1];
	  if (delta*delta > tolerance*tolerance)
	    return 0;
	}

      # check passed
      return 1;
    }
  else if (key == "CONTRIBUTION")
    {
      tolerance = 5e-4;

      # number of lines must match
      if (i1 != i2)
	return 0;

      # all pairs of lines must match
      for (j=0; j<i1; j++)
	{
	  # energy is the field after =
	  split(line1[j], x1);
	  split(line2[j], x2);
	  for (k in x1)
	    if (k > 2 && match(x1[k], "="))
	      break;
	  delta = x1[k+1] - x2[k+1];
	  if (delta*delta > tolerance*tolerance)
	    return 0;
	}

      # check passed
      return 1;
    }
  else if (key == "BANDS")
    {
      tolerance = 1e-1;

      # number of lines must match
      if (i1 != i2)
	return 0;

      # all pairs of lines must match
      for (j=0; j<i1; j++)
	{
	  # number of eigenvalues must match
	  f1 = split(line1[j], x1);
	  f2 = split(line2[j], x2);
	  if (f1 != f2)
	    return 0;

	  # all eigenvalues must match
	  for (k in x1)
	    if (k > 1)
	      {
		delta = x1[k] - x2[k];
		if (delta*delta > tolerance*tolerance)
		  return 0;
	      }
	}

      # check passed
      return 1;
    }
  else if (key == "PRESSURE")
    {
      tolerance = 1e-0;

      # there must be exactly one line
      if (i1 != 1 || i2 != 1)
	return 0;

      split(line1[0], x1);
      split(line2[0], x2);
      delta = x1[8] - x2[8];
      if (delta*delta > tolerance*tolerance)
	return 0;

      # check passed
      return 1;
    }
  else if (key == "STRESS")
    {
      tol_au   = 1e-5;
      tol_kbar = 1e-0;

      # number of lines must match
      if (i1 != i2)
	return 0;

      # all pairs of lines must match
      for (j=0; j<i1; j++)
	{
	  # all components must match
	  split(line1[j], x1);
	  split(line2[j], x2);
	  for (k=3; k<=5; k++)
	    {
	      delta = x1[k] - x2[k];
	      if (delta*delta > tol_au*tol_au)
		return 0;
	    }
	  for (k=6; k<=8; k++)
	    {
	      delta = x1[k] - x2[k];
	      if (delta*delta > tol_kbar*tol_kbar)
		return 0;
	    }
	}

      # check passed
      return 1;
    }
  else if (key == "DDV")
    {
      tolerance = 1e-3;

      # there must be exactly one line
      if (i1 != 1 || i2 != 1)
	return 0;

      split(line1[0], x1);
      split(line2[0], x2);
      delta = x1[10] - x2[10];
      if (delta*delta > tolerance*tolerance)
	return 0;

      # check passed
      return 1;
    }
  else if (key == "DIELECTRIC")
    {
      tolerance = 1e-6;

      # number of lines must match
      if (i1 != i2)
	return 0;

      # all pairs of lines must match
      for (j=0; j<i1; j++)
	{
	  # all components must match
	  split(line1[j], x1);
	  split(line2[j], x2);
	  for (k=4; k<=6; k++)
	    {
	      delta = x1[k] - x2[k];
	      if (delta*delta > tolerance*tolerance)
		return 0;
	    }
	}

      # check passed
      return 1;
    }
  else if (key == "OMEGA")
    {
      tol_thz = 1e-2;
      tol_cm  = 2e-1;

      # number of lines must match
      if (i1 != i2)
	return 0;

      # all pairs of lines must match
      for (j=0; j<i1; j++)
	{
	  split(line1[j], x1);
	  split(line2[j], x2);
	  for (k in x1)
	    if (k > 2 && match(x1[k], "THz"))
	      break;
	  delta1 = x1[k-1] - x2[k-1];
	  delta2 = x1[k+2] - x2[k+2];
	  if (delta1*delta1 > tol_thz*tol_thz \
	      || delta2*delta2 > tol_cm*tol_cm)
	    return 0;
	}

      # check passed
      return 1;
    }
  else if (key == "FORCE")
    {
      tolerance = 1e-4;

      # number of lines must match
      if (i1 != i2)
	return 0;

      # all pairs of lines must match
      for (j=0; j<i1; j++)
	{
	  # all components must match
	  split(line1[j], x1);
	  split(line2[j], x2);
	  for (k=9; k<=11; k++)
	    {
	      delta = x1[k] - x2[k];
	      if (delta*delta > tolerance*tolerance)
		return 0;
	    }
	}

      # check passed
      return 1;
    }
  else if (key == "TFORCE")
    {
      tol_force = 1e-4;
      tol_scf   = 1e-4;

      # there must be exactly one line
      if (i1 != 1 || i2 != 1)
	return 0;

      split(line1[0], x1);
      split(line2[0], x2);
      delta1 = x1[6] - x2[6];
      delta2 = x1[11] - x2[11];
      if (delta1*delta1 > tol_force*tol_force \
	  || delta2*delta2 > tol_scf*tol_scf)
	return 0;

      # check passed
      return 1;
    }
  else if (key == "TEMPERATURE")
    {
      tol_ekin = 1e-6;
      tol_t    = 2e-1;
      tol_etot = 1e-5;

      # there must be exactly one line
      if (i1 != 1 || i2 != 1)
	return 0;

      split(line1[0], x1);
      split(line2[0], x2);
      delta1 = x1[5] - x2[5];
      delta2 = x1[9] - x2[9];
      delta3 = x1[13] - x2[13];
      if (delta1*delta1 > tol_ekin*tol_ekin || delta2*delta2 > tol_t*tol_t \
	  || delta3*delta3 > tol_etot*tol_etot)
	return 0;

      # check passed
      return 1;
    }
  else if (key == "POSITIONS")
    {
      tolerance = 1e-3;

      # number of lines must match
      if (i1 != i2)
	return 0;

      # all pairs of lines must match
      for (j=0; j<i1; j++)
	{
	  # all components must match
	  split(line1[j], x1);
	  split(line2[j], x2);
	  for (k=4; k<=6; k++)
	    {
	      delta = x1[k] - x2[k];
	      if (delta*delta > tolerance*tolerance)
		return 0;
	    }
	}

      # check passed
      return 1;
    }
  else if (key == "EFINAL")
    {
      tolerance = 1e-6;

      # there must be exactly one line
      if (i1 != 1 || i2 != 1)
	return 0;

      split(line1[0], x1);
      split(line2[0], x2);
      delta = x1[7] - x2[7];
      if (delta*delta > tolerance*tolerance)
	return 0;

      # check passed
      return 1;
    }
  else if (key == "CHECKPOINT")
    {
      return 1;
    }
  else
    {
      # no key

      # number of lines must match
      if (i1 != i2)
	return 0;

      # all pairs of lines must match
      for (j=0; j<i1; j++)
	{
	  split(line1[j], x1);
	  split(line2[j], x2);
	  for (k in x1)
	    if (k > 1 && x1[k] != x2[k])
	      return 0;
	}

      # check passed
      return 1;
    }
}

function print_diff()
{
  # recreate output of diff
  print_head();
  for (i=0; i<i1; i++)
    print strip_key(line1[i]);
  if (i1 > 0 && i2 > 0)
    print "---";
  for (i=0; i<i2; i++)
    print strip_key(line2[i]);
}

function read_head()
{
  # read a header line from output of diff
  # s1 is the first non-matching line in the first file
  # n1 is the number of non-matching lines
  # s2, n2 are the same things for the second file

  type = $0; gsub("[0-9,]", "", type);
  split($0, x, type);

  if (type == "a")
    {
      s1 = x[1]; n1 = 0;
    }
  else if (match(x[1], ","))
    {
      split(x[1], y, ",");
      s1 = y[1]; n1 = y[2] - s1 + 1;
    }
  else
    {
      s1 = x[1]; n1 = 1;
    }

  if (type == "d")
    {
      s2 = x[2]; n2 = 0;
    }
  else if (match(x[2], ","))
    {
      split(x[2], y, ",");
      s2 = y[1]; n2 = y[2] - s2 + 1;
    }
  else
    {
      s2 = x[2]; n2 = 1;
    }
}

function print_head()
{
  # o1, o2 must be given as command line arguments

  if (i1 == 0)
    {
      if (i2 == 0)
	; # nothing
      else if (i2 == 1)
	printf("%da%d\n", s1+o1, s2+o2);
      else
	printf("%da%d,%d\n", s1+o1, s2+o2, s2+o2+i2-1);
    }
  else if (i1 == 1)
    {
      if (i2 == 0)
	printf("%dd%d\n", s1+o1, s2+o2);
      else if (i2 == 1)
	printf("%dc%d\n", s1+o1, s2+o2);
      else
	printf("%dc%d,%d\n", s1+o1, s2+o2, s2+o2+i2-1);
    }
  else
    {
      if (i2 == 0)
	printf("%d,%dd%d\n", s1+o1, s1+o1+i1-1, s2+o2);
      else if (i2 == 1)
	printf("%d,%dc%d\n", s1+o1, s1+o1+i1-1, s2+o2);
      else
	printf("%d,%dc%d,%d\n", s1+o1, s1+o1+i1-1, s2+o2, s2+o2+i2-1);
    }
}

function get_key(line)
{
  if (match(line, "@.*@"))
    {
      split(line, x, "@");
      return x[2];
    }
  else
    return "";
}

function strip_key(line)
{
  sub("@.*@ ", "", line);
  return line;
}
