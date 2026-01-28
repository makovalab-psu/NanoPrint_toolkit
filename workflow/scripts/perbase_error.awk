# Script.awk
BEGIN {
    # Create a lookup array for ASCII values of characters
    for (n = 0; n < 256; n++) {
        ord[sprintf("%c", n)] = n;  # Map each character to its ASCII value
    }
}
{
    # Initialize variables
    sum_error_prob = 0;
    count = 1;

    # Split the quality score string (column 6) into individual characters
    split($6, qualities, "");

    # Loop through each character in the quality score string
    for (i = 2; i <= length($6); i++) {  # Start from 2 to skip the first empty element
        # Get the character from the quality score string
        char = substr($6, i, 1);
        
        # Convert character to ASCII value and then to Phred score
        ascii_value = ord[char];
        phred_score = ascii_value - 33; 
	

        # Calculate error probability
        error_prob = 10^(-phred_score / 10);

        # Accumulate the error probability and count
        sum_error_prob += error_prob;
        count++;
    }

    # Calculate mean error probability
    mean_error_prob = (count > 0) ? sum_error_prob / count : 0;

    # Print the original line with the mean error probability
    print $1 "\t" $2 "\t" $3 "\t" count "\t"  mean_error_prob
}
